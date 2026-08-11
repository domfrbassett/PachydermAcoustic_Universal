using HDF.PInvoke;
using System;
using System.Linq;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using System.Text;
using Hare.Geometry;
using Vector = Hare.Geometry.Vector;
using Pachyderm_Acoustic.Utilities;
using Pachyderm_Acoustic.Audio;

namespace Pachyderm_Acoustic
{
    namespace Audio
    {
        public class HRTF
        {
            public Topology T;
            public Voxel_Grid VG;
            List<Vector> emitterVectors;
            double[,] emitterdirs;
            float[][][] HRIR;
            private int Fs;
            Vector[] PrincipalDirections;

            public bool AngularDistancePass { get; private set; }
            public bool GlobeCoveragePass { get; private set; }
            public string ValidationMessage { get; private set; }
            public bool SofaCompatibilityPass { get; private set; }
            public bool ValidationPassed => SofaCompatibilityPass;

            public double AvgAngularDistanceFrontal { get; private set; }
            public double AvgAngularDistanceRear { get; private set; }
            public double MaxCoverageGap { get; private set; }

            double[][][] Loaded_HRIR;
            double[][][] Loaded_HRIR_44100;

            public HRTF(string SofafilePath)
            {
                // Open the .SOFA file
                long fileId = H5F.open(SofafilePath, H5F.ACC_RDONLY);
                if (fileId < 0)
                {
                    throw new InvalidOperationException("Error opening SOFA file.");
                }

                try
                {
                    string Conv = ReadGlobalString(fileId, "Conventions");
                    string Conv_SOFA = ReadGlobalString(fileId, "SOFAConventions");
                    string DataType = ReadGlobalString(fileId, "DataType");
                    ValidateSofaCompatibility(Conv, Conv_SOFA, DataType);

                    emitterdirs = ReadSourcePosition(fileId);
                    if (emitterdirs == null || emitterdirs.GetLength(0) == 0)
                        throw new InvalidOperationException("No usable SourcePosition directions were found in the SOFA file.");

                    Fs = ReadSamplingFrequency(fileId);
                    ReadHrtfDataset(fileId, "Data.IR", emitterdirs.GetLength(0));
                    BuildDirectionalSamplingSummary();

                    emitterVectors = new List<Vector>();
                    for (int i = 0; i < emitterdirs.GetLength(0); i++)
                    {
                        double azRad = Math.PI * emitterdirs[i, 0] / 180.0;
                        double elRad = Math.PI * emitterdirs[i, 1] / 180.0;
                        double x = Math.Cos(elRad) * Math.Cos(azRad);
                        double y = Math.Cos(elRad) * Math.Sin(azRad);
                        double z = Math.Sin(elRad);
                        Vector d = new Vector(x, y, z);
                        d.Normalize();
                        emitterVectors.Add(d);
                    }

                    var validHRIRs = new List<double[][]>();
                    var validDirections = new List<Vector>();

                    for (int i = 0; i < emitterVectors.Count; i++)
                    {
                        validHRIRs.Add(Pach_SP_HRTF.ResampleHRIRWDL(HRIR[i], Fs, 44100));
                        validDirections.Add(emitterVectors[i]);
                    }

                    Loaded_HRIR_44100 = validHRIRs.ToArray(); // Store HRIRs resampled to 44100 Hz - useful for many applications and can be retrieved quickly.
                    Directions = validDirections.ToArray();
                }
                finally
                {
                    H5F.close(fileId);
                }
            }

            private void ValidateSofaCompatibility(string conventions, string sofaConvention, string dataType)
            {
                if (!string.Equals((conventions ?? string.Empty).Trim(), "SOFA", StringComparison.OrdinalIgnoreCase))
                    throw new InvalidOperationException("The selected file does not declare the SOFA convention.");

                string conv = (sofaConvention ?? string.Empty).Trim();
                string type = (dataType ?? string.Empty).Trim();

                if (!string.IsNullOrWhiteSpace(type) && !string.Equals(type, "FIR", StringComparison.OrdinalIgnoreCase))
                    throw new InvalidOperationException($"Unsupported SOFA data type '{type}'. Pachyderm currently expects FIR HRIR data stored in Data.IR.");

                if (string.Equals(conv, "SimpleFreeFieldHRIR", StringComparison.OrdinalIgnoreCase))
                {
                    SofaCompatibilityPass = true;
                    return;
                }

                if (string.Equals(conv, "GeneralFIR", StringComparison.OrdinalIgnoreCase))
                {
                    SofaCompatibilityPass = true;
                    ValidationMessage = "Warning: this file uses GeneralFIR. Pachyderm will treat SourcePosition as HRTF directions and Data.IR channels 0 and 1 as ears.";
                    return;
                }

                if (string.Equals(conv, "SimpleFreeFieldHRTF", StringComparison.OrdinalIgnoreCase) ||
                    string.Equals(conv, "FreeFieldHRTF", StringComparison.OrdinalIgnoreCase))
                    throw new InvalidOperationException($"Unsupported SOFA convention '{conv}'. Frequency-domain HRTF conventions require Data.Real/Data.Imag conversion before Pachyderm can use them.");

                if (string.Equals(conv, "SimpleFreeFieldHRSOS", StringComparison.OrdinalIgnoreCase) ||
                    string.Equals(conv, "SimpleFreeFieldSOS", StringComparison.OrdinalIgnoreCase))
                    throw new InvalidOperationException($"Unsupported SOFA convention '{conv}'. SOS HRTF conventions must be converted to FIR HRIR data before Pachyderm can use them.");

                if (string.Equals(conv, "FreeFieldHRIR", StringComparison.OrdinalIgnoreCase) ||
                    string.Equals(conv, "GeneralFIR-E", StringComparison.OrdinalIgnoreCase))
                    throw new InvalidOperationException($"Unsupported SOFA convention '{conv}'. Emitter-dependent FIR data is not handled by this HRTF loader.");

                throw new InvalidOperationException($"Unsupported SOFA convention '{(string.IsNullOrWhiteSpace(conv) ? "unknown" : conv)}'. Pachyderm currently supports SimpleFreeFieldHRIR and compatible GeneralFIR files.");
            }

            private void BuildDirectionalSamplingSummary()
            {
                var vectors = Enumerable.Range(0, emitterdirs.GetLength(0))
                    .Select(i => PachTools.SphericalToCartesian(emitterdirs[i, 0], emitterdirs[i, 1]))
                    .ToArray();

                var frontVectors = vectors.Where(v => v.dx >= 0).ToArray();
                var backVectors = vectors.Where(v => v.dx < 0).ToArray();
                AngularDistancePass = frontVectors.Length > 0 && backVectors.Length > 0;
                AvgAngularDistanceFrontal = frontVectors.Length > 1 ? AverageMinAngularDistance(frontVectors) : 180.0;
                AvgAngularDistanceRear = backVectors.Length > 1 ? AverageMinAngularDistance(backVectors) : 180.0;

                MaxCoverageGap = EstimateMaxCoverageGap(vectors, 362);
                bool hasUpper = emitterdirs.Cast<double>().Where((_, i) => i % 2 == 1).Any(el => el > 35.0);
                bool hasLower = emitterdirs.Cast<double>().Where((_, i) => i % 2 == 1).Any(el => el < -35.0);
                GlobeCoveragePass = hasUpper && hasLower && MaxCoverageGap <= 60.0;

                var warnings = new List<string>();
                if (!string.IsNullOrWhiteSpace(ValidationMessage)) warnings.Add(ValidationMessage);
                if (vectors.Length < 24) warnings.Add("Warning: the SOFA file has very few HRTF directions; binaural rendering may be coarse.");
                if (!AngularDistancePass) warnings.Add("Warning: SourcePosition does not cover both front and rear hemispheres.");
                if (!hasUpper || !hasLower) warnings.Add("Warning: SourcePosition does not appear to cover both upper and lower hemispheres.");
                if (MaxCoverageGap > 60.0) warnings.Add($"Warning: estimated maximum directional gap is {MaxCoverageGap:F1} deg; full-3D binaural rendering may be sparse.");

                ValidationMessage = warnings.Count > 0 ? string.Join(System.Environment.NewLine, warnings) : string.Empty;
            }

            private static double EstimateMaxCoverageGap(Vector[] directions, int sampleCount)
            {
                if (directions == null || directions.Length == 0) return 180.0;

                double maxGap = 0.0;
                double goldenAngle = Math.PI * (3.0 - Math.Sqrt(5.0));

                for (int i = 0; i < sampleCount; i++)
                {
                    double z = 1.0 - 2.0 * (i + 0.5) / sampleCount;
                    double r = Math.Sqrt(Math.Max(0.0, 1.0 - z * z));
                    double theta = i * goldenAngle;
                    var sample = new Vector(Math.Cos(theta) * r, Math.Sin(theta) * r, z);

                    double minAngle = directions.Min(direction => AngularDistanceDegrees(sample, direction));
                    if (minAngle > maxGap) maxGap = minAngle;
                }

                return maxGap;
            }

            public static int GetRequiredSubdivisionOrder(int totalSources)
            {
                if (totalSources <= 20) return 0;
                if (totalSources <= 80) return 1;
                if (totalSources <= 320) return 2;
                return 3;
            }
            public int SampleCt
            {
                get
                {
                    return HRIR[0][0].Length;
                }
            }

            public int DirsCT
            {
                get
                {
                    return Directions == null ? 0 : Directions.Length;
                }
            }

            public Vector[] Directions
            {
                get { return PrincipalDirections; }
                set { PrincipalDirections = value; }
            }

            private bool CheckHrtfAngularDistance()
            {
                int nSources = emitterdirs.GetLength(0);

                // Convert spherical to cartesian unit vectors
                var vectors = new Vector[nSources];
                for (int i = 0; i < nSources; i++)
                {
                    vectors[i] = PachTools.SphericalToCartesian(emitterdirs[i, 0], emitterdirs[i, 1]);
                }

                // Separate front/back using X coordinate (positive = front, negative = back)
                var frontVectors = vectors.Where(v => v.dx >= 0).ToArray();
                var backVectors = vectors.Where(v => v.dx < 0).ToArray();

                // Fail if either hemisphere has no sources
                if (frontVectors.Length == 0 || backVectors.Length == 0)
                    return false;

                // Compute average minimum angular distance
                AvgAngularDistanceFrontal = AverageMinAngularDistance(frontVectors);
                AvgAngularDistanceRear = AverageMinAngularDistance(backVectors);

                // Check against threshold (8°)
                return AvgAngularDistanceFrontal <= 8 && AvgAngularDistanceRear <= 8;
            }

            private bool CheckHrtfGlobeCoverage(double maxAllowedGapDeg)
            {
                var emitters = Enumerable.Range(0, emitterdirs.GetLength(0))
                    .Select(i => PachTools.SphericalToCartesian(emitterdirs[i, 0], emitterdirs[i, 1]))
                    .ToArray();

                int N_az = 73; // Number of azimuthal divisions - approximately 5° spacing
                int N_el = 19; // Number of elevation divisions - approximately 5° spacing

                var refPoints = new List<Hare.Geometry.Vector>();

                for (int elIdx = 0; elIdx < N_el; elIdx++)
                {
                    double el = -90.0 + 180.0 * elIdx / (N_el - 1);
                    for (int azIdx = 0; azIdx < N_az; azIdx++)
                    {
                        double az = 360.0 * azIdx / N_az;
                        refPoints.Add(PachTools.SphericalToCartesian(az, el));
                    }
                }

                MaxCoverageGap = 0.0;
                for (int i = 0; i < refPoints.Count; i++)
                {
                    double minAngle = emitters.Min(emitter => AngularDistanceDegrees(refPoints[i], emitter));
                    if (minAngle > MaxCoverageGap) MaxCoverageGap = minAngle;

                    if (minAngle > maxAllowedGapDeg)
                    {
                        return false;
                    }
                }

                return true;
            }

            public static double AngularDistanceDegrees(Hare.Geometry.Vector v1, Hare.Geometry.Vector v2)
            {
                double dot = Hare.Geometry.Hare_math.Dot(v1, v2);
                dot = Math.Min(1.0, Math.Max(-1.0, dot)); // Clamp to [-1, 1] to avoid NaN
                return Math.Acos(dot) * 180.0 / Math.PI; // Convert radians to degrees
            }

            private double AverageMinAngularDistance(Hare.Geometry.Vector[] vectors)
            {
                int n = vectors.Length;
                double[] minAngles = new double[n];

                for (int i = 0; i < n; i++)
                {
                    double minAngle = 180.0;
                    for (int j = 0; j < n; j++)
                    {
                        if (i == j) continue;
                        double angle = AngularDistanceDegrees(vectors[i], vectors[j]);
                        if (angle < minAngle) minAngle = angle;
                    }
                    minAngles[i] = minAngle;
                }

                return minAngles.Average();
            }

            public double GetAvgAngularDistanceFrontal()
            {
                return AvgAngularDistanceFrontal;
            }

            public double GetAvgAngularDistanceRear()
            {
                return AvgAngularDistanceRear;
            }

            public double GetMaxCoverageGap()
            {
                return MaxCoverageGap;
            }

            public int ReadSamplingFrequency(long fileId)
            {
                long datasetId = H5D.open(fileId, "/Data.SamplingRate");
                if (datasetId < 0)
                {
                    throw new Exception("Failed to open dataset '/Data.SamplingRate'. This dataset is required to get the sampling frequency.");
                }

                try
                {
                    long typeId = H5D.get_type(datasetId);
                    if (typeId < 0)
                    {
                        throw new Exception("Failed to get datatype for '/Data.SamplingRate'.");
                    }

                    try
                    {
                        double[] Fs = new double[1];
                        GCHandle hnd = GCHandle.Alloc(Fs, GCHandleType.Pinned);

                        try
                        {
                            IntPtr ptr = hnd.AddrOfPinnedObject();
                            if (H5D.read(datasetId, typeId, H5S.ALL, H5S.ALL, H5P.DEFAULT, ptr) < 0)
                            {
                                throw new Exception("Failed to read '/Data.SamplingRate'.");
                            }

                            System.Diagnostics.Debug.WriteLine($"Input sampling rate: {(int)Math.Round(Fs[0])} Hz");
                            System.Diagnostics.Debug.WriteLine("Output sampling rate (target): 44100 Hz");
                            return (int)Math.Round(Fs[0]);
                        }
                        finally
                        {
                            hnd.Free();
                        }
                    }
                    finally
                    {
                        H5T.close(typeId);
                    }
                }
                finally
                {
                    H5D.close(datasetId);
                }
            }

            private double[,] ReadSourcePosition(long fileId)
            {
                string[] candidates = { "SourcePosition" };

                foreach (var name in candidates)
                {
                    long datasetId = H5D.open(fileId, name);
                    if (datasetId < 0)
                    {
                        Console.WriteLine($"Dataset {name} not found.");
                        continue;
                    }

                    long spaceId = -1;

                    try
                    {
                        spaceId = H5D.get_space(datasetId);
                        if (spaceId < 0)
                        {
                            Console.WriteLine($"Error getting dataspace for {name}.");
                            continue;
                        }

                        ulong[] dims = new ulong[2];
                        int rank = H5S.get_simple_extent_dims(spaceId, dims, null);
                        if (rank != 2)
                        {
                            Console.WriteLine($"Unexpected {name} rank. Expected 2D Nx3 or 3xN, got rank {rank}.");
                            continue;
                        }

                        bool transposed = false;
                        int numSources;
                        if (dims[1] == 3)
                        {
                            numSources = (int)dims[0];
                        }
                        else if (dims[0] == 3)
                        {
                            transposed = true;
                            numSources = (int)dims[1];
                        }
                        else
                        {
                            Console.WriteLine($"Unexpected {name} format. Expected Nx3 or 3xN, got {dims[0]}x{dims[1]}.");
                            continue;
                        }

                        int len = numSources * 3;
                        double[] flat = new double[len];

                        GCHandle h = GCHandle.Alloc(flat, GCHandleType.Pinned);
                        try
                        {
                            IntPtr ptr = h.AddrOfPinnedObject();

                            if (H5D.read(datasetId, H5T.NATIVE_DOUBLE, H5S.ALL, H5S.ALL, H5P.DEFAULT, ptr) < 0)
                            {
                                Console.WriteLine($"Failed to read {name} dataset.");
                                continue;
                            }
                        }
                        finally
                        {
                            h.Free();
                        }

                        string positionType = ReadDatasetStringAttribute(fileId, name, "Type");
                        string units = ReadDatasetStringAttribute(fileId, name, "Units");
                        bool isCartesian = !string.IsNullOrWhiteSpace(positionType) && positionType.IndexOf("cartesian", StringComparison.OrdinalIgnoreCase) >= 0;
                        bool isSpherical = string.IsNullOrWhiteSpace(positionType) || positionType.IndexOf("spherical", StringComparison.OrdinalIgnoreCase) >= 0;

                        if (!isSpherical && !isCartesian)
                        {
                            Console.WriteLine($"Unsupported {name}:Type '{positionType}'. Expected spherical or cartesian.");
                            continue;
                        }

                        bool angularUnitsAreRadians = !string.IsNullOrWhiteSpace(units) && units.IndexOf("radian", StringComparison.OrdinalIgnoreCase) >= 0;
                        double[,] result = new double[numSources, 2];

                        for (int i = 0; i < numSources; ++i)
                        {
                            double v0 = transposed ? flat[i] : flat[i * 3 + 0];
                            double v1 = transposed ? flat[numSources + i] : flat[i * 3 + 1];
                            double v2 = transposed ? flat[2 * numSources + i] : flat[i * 3 + 2];

                            if (isCartesian)
                            {
                                double xy = Math.Sqrt(v0 * v0 + v1 * v1);
                                double r = Math.Sqrt(v0 * v0 + v1 * v1 + v2 * v2);
                                if (r <= 1e-12)
                                    throw new InvalidOperationException($"{name} contains a zero-length cartesian source vector at index {i}.");

                                result[i, 0] = Math.Atan2(v1, v0) * 180.0 / Math.PI;
                                result[i, 1] = Math.Atan2(v2, xy) * 180.0 / Math.PI;
                            }
                            else
                            {
                                result[i, 0] = angularUnitsAreRadians ? v0 * 180.0 / Math.PI : v0;
                                result[i, 1] = angularUnitsAreRadians ? v1 * 180.0 / Math.PI : v1;
                            }
                        }

                        Console.WriteLine($"Read {numSources} source positions from {name} ({(isCartesian ? "cartesian" : "spherical")}).");
                        return result;
                    }
                    finally
                    {
                        if (spaceId >= 0) H5S.close(spaceId);
                        if (datasetId >= 0) H5D.close(datasetId);
                    }
                }

                Console.WriteLine("No valid SourcePosition dataset found.");
                return null;
            }

            private string ReadDatasetStringAttribute(long fileId, string datasetName, string attributeName)
            {
                long attrId = H5A.open_by_name(fileId, datasetName, attributeName, H5P.DEFAULT, H5P.DEFAULT);
                if (attrId < 0) return null;

                try
                {
                    long typeId = H5A.get_type(attrId);
                    if (typeId < 0) return null;

                    try
                    {
                        int size = Math.Max(H5T.get_size(typeId).ToInt32(), 1);
                        byte[] buffer = new byte[size];
                        GCHandle pinnedArray = GCHandle.Alloc(buffer, GCHandleType.Pinned);
                        try
                        {
                            if (H5A.read(attrId, typeId, pinnedArray.AddrOfPinnedObject()) < 0)
                                return null;
                        }
                        finally
                        {
                            pinnedArray.Free();
                        }

                        return Encoding.ASCII.GetString(buffer).TrimEnd('\0', ' ');
                    }
                    finally
                    {
                        H5T.close(typeId);
                    }
                }
                finally
                {
                    H5A.close(attrId);
                }
            }
            private string ReadGlobalString(long fileId, string Field)
            {
                // Open the attribute
                long attrId = H5A.open_by_name(fileId, "/", Field, H5P.DEFAULT, H5P.DEFAULT);
                if (attrId < 0)
                {
                    return null;
                }

                // Get the datatype of the attribute
                long typeId = H5A.get_type(attrId);
                if (typeId < 0)
                {
                    H5A.close(attrId);
                    return null;
                }

                // Determine the size of the attribute
                IntPtr size = H5T.get_size(typeId);
                byte[] buffer = new byte[size.ToInt32()];

                // Read the attribute
                GCHandle pinnedArray = GCHandle.Alloc(buffer, GCHandleType.Pinned);
                IntPtr pointer = pinnedArray.AddrOfPinnedObject();
                if (H5A.read(attrId, typeId, pointer) < 0)
                {
                    Console.WriteLine("Error reading attribute.");
                    pinnedArray.Free();
                    H5T.close(typeId);
                    H5A.close(attrId);
                    return null;
                }
                pinnedArray.Free();

                // Convert the byte array to a string
                string conventions = Encoding.ASCII.GetString(buffer).TrimEnd('\0');

                // Close resources
                H5T.close(typeId);
                H5A.close(attrId);

                return conventions;
            }

            private float[] ReadSourcePositionUnits(long fileId)
            {
                // Open the attribute
                long attrId = H5A.open_by_name(fileId, "SourcePosition", "Units", H5P.DEFAULT, H5P.DEFAULT);
                if (attrId < 0)
                {
                    return null;
                }

                // Get the datatype of the attribute
                long typeId = H5A.get_type(attrId);
                if (typeId < 0)
                {
                    H5A.close(attrId);
                    return null;
                }

                // Get the dataspace of the attribute
                long spaceId = H5A.get_space(attrId);
                if (spaceId < 0)
                {
                    H5T.close(typeId);
                    H5A.close(attrId);
                    return null;
                }

                // Determine the number of elements in the dataspace
                ulong[] dims = new ulong[1];
                int rank = H5S.get_simple_extent_dims(spaceId, dims, null);
                if (rank < 0)
                {
                    Console.WriteLine("Error getting dataspace dimensions.");
                    H5S.close(spaceId);
                    H5T.close(typeId);
                    H5A.close(attrId);
                    return null;
                }

                // Allocate a buffer for the float array
                float[] buffer = new float[dims[0]];

                // Read the attribute
                GCHandle pinnedArray = GCHandle.Alloc(buffer, GCHandleType.Pinned);
                IntPtr pointer = pinnedArray.AddrOfPinnedObject();
                if (H5A.read(attrId, typeId, pointer) < 0)
                {
                    Console.WriteLine("Error reading attribute.");
                    pinnedArray.Free();
                    H5S.close(spaceId);
                    H5T.close(typeId);
                    H5A.close(attrId);
                    return null;
                }
                pinnedArray.Free();

                // Close resources
                H5S.close(spaceId);
                H5T.close(typeId);
                H5A.close(attrId);

                return buffer;
            }

            private void ReadHrtfDataset(long fileId, string datasetName, int expectedMeasurements)
            {
                long datasetId = H5D.open(fileId, datasetName);
                if (datasetId < 0)
                    throw new InvalidOperationException($"The SOFA file does not contain required FIR dataset {datasetName}.");

                try
                {
                    long dataspaceId = H5D.get_space(datasetId);
                    try
                    {
                        int rank = H5S.get_simple_extent_ndims(dataspaceId);
                        if (rank != 3)
                            throw new InvalidOperationException($"Unsupported {datasetName} rank {rank}. Pachyderm expects SimpleFreeFieldHRIR-style M x R x N FIR data.");

                        ulong[] dims = new ulong[3];
                        if (H5S.get_simple_extent_dims(dataspaceId, dims, null) < 0)
                            throw new InvalidOperationException($"Unable to read dimensions for {datasetName}.");

                        int measurements = (int)dims[0];
                        int channels = (int)dims[1];
                        int samples = (int)dims[2];

                        if (measurements != expectedMeasurements)
                            throw new InvalidOperationException($"SourcePosition count ({expectedMeasurements}) does not match {datasetName} measurement count ({measurements}).");
                        if (channels < 2)
                            throw new InvalidOperationException($"{datasetName} has {channels} receiver channel(s); binaural rendering requires at least two.");
                        if (samples <= 0)
                            throw new InvalidOperationException($"{datasetName} contains no FIR samples.");

                        float[,,] dataIR = new float[measurements, channels, samples];
                        GCHandle handle = GCHandle.Alloc(dataIR, GCHandleType.Pinned);
                        try
                        {
                            if (H5D.read(datasetId, H5T.NATIVE_FLOAT, H5S.ALL, H5S.ALL, H5P.DEFAULT, handle.AddrOfPinnedObject()) < 0)
                                throw new InvalidOperationException($"Failed to read {datasetName}.");
                        }
                        finally
                        {
                            handle.Free();
                        }

                        HRIR = new float[measurements][][];
                        for (int i = 0; i < measurements; i++)
                        {
                            HRIR[i] = new float[channels][];
                            for (int ch = 0; ch < channels; ch++)
                            {
                                HRIR[i][ch] = new float[samples];
                                for (int k = 0; k < samples; k++)
                                    HRIR[i][ch][k] = dataIR[i, ch, k];
                            }
                        }
                    }
                    finally
                    {
                        H5S.close(dataspaceId);
                    }
                }
                finally
                {
                    H5D.close(datasetId);
                }
            }
            double[][] Loaded_Filter;

            public void Load(Direct_Sound Direct, ImageSourceData ISData, Pachyderm_Acoustic.Environment.Receiver_Bank RTData, SystemResponseCompensation.SystemCompensationSettings sysCompSettings, double CO_Time_ms, int targetFs, int Rec_ID, bool Start_at_Zero, bool flat, bool auto)
            {
                Loaded_Filter = new double[6][];
                Loaded_Filter[0] = Pachyderm_Acoustic.Utilities.IR_Construction.Aurfilter_Directional(Direct, ISData, RTData, CO_Time_ms, targetFs, Rec_ID, Start_at_Zero, 0, 0, true, flat);
                Loaded_Filter[1] = Pachyderm_Acoustic.Utilities.IR_Construction.Aurfilter_Directional(Direct, ISData, RTData, CO_Time_ms, targetFs, Rec_ID, Start_at_Zero, 0, 180, true, flat);
                Loaded_Filter[2] = Pachyderm_Acoustic.Utilities.IR_Construction.Aurfilter_Directional(Direct, ISData, RTData, CO_Time_ms, targetFs, Rec_ID, Start_at_Zero, 0, 90, true, flat);
                Loaded_Filter[3] = Pachyderm_Acoustic.Utilities.IR_Construction.Aurfilter_Directional(Direct, ISData, RTData, CO_Time_ms, targetFs, Rec_ID, Start_at_Zero, 0, 270, true, flat);
                Loaded_Filter[4] = Pachyderm_Acoustic.Utilities.IR_Construction.Aurfilter_Directional(Direct, ISData, RTData, CO_Time_ms, targetFs, Rec_ID, Start_at_Zero, 90, 0, true, flat);
                Loaded_Filter[5] = Pachyderm_Acoustic.Utilities.IR_Construction.Aurfilter_Directional(Direct, ISData, RTData, CO_Time_ms, targetFs, Rec_ID, Start_at_Zero, -90, 0, true, flat);

                if (targetFs != Fs)
                {
                    if (targetFs == 44100)
                    {
                        var validHRIRs = new List<double[][]>();

                        for (int i = 0; i < Loaded_HRIR_44100.Length; i++)
                        {
                            var hrirCopy = new double[Loaded_HRIR_44100[i].Length][];
                            for (int ch = 0; ch < Loaded_HRIR_44100[i].Length; ch++)
                            {
                                hrirCopy[ch] = new double[Loaded_HRIR_44100[i][ch].Length];
                                Array.Copy(Loaded_HRIR_44100[i][ch], hrirCopy[ch], Loaded_HRIR_44100[i][ch].Length);
                            }
                            validHRIRs.Add(hrirCopy);
                        }

                        Loaded_HRIR = validHRIRs.ToArray();
                    }
                    else
                    {
                        var validHRIRs = new List<double[][]>();
                        for (int i = 0; i < DirsCT; i++)
                        {
                            validHRIRs.Add(Pach_SP_HRTF.ResampleHRIRWDL(
                                HRIR[i],
                                Fs,
                                targetFs
                            ));
                        }
                        Loaded_HRIR = validHRIRs.ToArray();
                    }
                }
                else
                {
                    var validHRIRs = new List<double[][]>();

                    for (int i = 0; i < DirsCT; i++)
                    {
                        var hrirCopy = new double[HRIR[i].Length][];

                        for (int ch = 0; ch < HRIR[i].Length; ch++)
                        {
                            hrirCopy[ch] = new double[HRIR[i][ch].Length];
                            for (int k = 0; k < HRIR[i][ch].Length; k++)
                            {
                                hrirCopy[ch][k] = HRIR[i][ch][k];
                            }
                        }

                        validHRIRs.Add(hrirCopy);
                    }

                    Loaded_HRIR = validHRIRs.ToArray();
                }

                Pach_SP_HRTF.ApplySystemCompensation(Loaded_HRIR, Directions, targetFs, 0, sysCompSettings, auto);
            }
            public double[][] Binaural_IR(double _azi, double _alt)
            {
                if (Loaded_Filter == null || Loaded_Filter.Length < 6)
                    throw new InvalidOperationException("HRTF filters have not been loaded.");

                int signalLength = Loaded_Filter.First(f => f != null).Length;
                double[][] directionalSignals = new double[Directions.Length][];
                for (int i = 0; i < directionalSignals.Length; i++)
                    directionalSignals[i] = new double[signalLength];

                List<int> activeDirections = new List<int>();

                Vector[] worldAxes = new Vector[]
                {
                    new Vector(1, 0, 0),
                    new Vector(-1, 0, 0),
                    new Vector(0, 1, 0),
                    new Vector(0, -1, 0),
                    new Vector(0, 0, 1),
                    new Vector(0, 0, -1)
                };

                for (int axis = 0; axis < worldAxes.Length; axis++)
                {
                    Vector headDirection = PachTools.Rotate_Vector(worldAxes[axis], _azi, _alt, true);
                    headDirection.Normalize();

                    List<Pach_SP_HRTF.DirectionalGain> gains = Pach_SP_HRTF.VbapGains(headDirection, Directions);
                    foreach (Pach_SP_HRTF.DirectionalGain gain in gains)
                    {
                        AddScaled(directionalSignals[gain.Index], Loaded_Filter[axis], gain.Gain);
                        if (!activeDirections.Contains(gain.Index)) activeDirections.Add(gain.Index);
                    }
                }

                if (activeDirections.Count == 0)
                    throw new InvalidOperationException("No active HRTF directions were selected for binaural rendering.");

                double[][] activeDirectionalSignals = activeDirections.Select(idx => directionalSignals[idx]).ToArray();
                double[] dryDirectionalSignal = Pach_SP_HRTF.BuildDrySignal(activeDirectionalSignals);
                double dryRMS = Pach_SP_HRTF.ComputeRMS(dryDirectionalSignal);

                double[][][] convolvedSignals = new double[activeDirections.Count][][];
                for (int i = 0; i < activeDirections.Count; i++)
                {
                    int idx = activeDirections[i];
                    convolvedSignals[i] = new double[2][];
                    convolvedSignals[i][0] = Pach_SP.FFT_Convolution_double(directionalSignals[idx], Loaded_HRIR[idx][0], 0);
                    convolvedSignals[i][1] = Pach_SP.FFT_Convolution_double(directionalSignals[idx], Loaded_HRIR[idx][1], 0);
                }

                double[][] Signal = Pach_SP_HRTF.SumAcrossDirections(convolvedSignals);
                Pach_SP_HRTF.NormaliseStereoByDryRMS(Signal, dryRMS);

                return Signal;
            }
        
            private static void AddScaled(double[] destination, double[] source, double gain)
            {
                if (destination == null || source == null) return;
                int length = Math.Min(destination.Length, source.Length);
                for (int i = 0; i < length; i++)
                    destination[i] += source[i] * gain;
            }
        }
    }
}
