package au.edu.adelaide.fxmr.model;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

import au.edu.adelaide.fxmr.data.GeneralModel;
import au.edu.adelaide.fxmr.model.mr.MRUtil;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;

public class ChealTests {
	private static final double NaN = Double.NaN;
	private double[][][] data = {
			{ { 0.58, 33.77, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN }, { -1.98, 17.47, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN },
					{ -0.36, 15.55, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN }, { 0.34, 16.28, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN },
					{ 0.7, 22.81, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN }, { 0.21, 21.5, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN },
					{ -0.89, 18.44, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN }, { 4.02, 14.64, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN },
					{ NaN, NaN, 2.32, 20.05, NaN, NaN, NaN, NaN, NaN, NaN }, { NaN, NaN, -1.39, 9.19, NaN, NaN, NaN, NaN, NaN, NaN },
					{ NaN, NaN, -1, 13.7, NaN, NaN, NaN, NaN, NaN, NaN }, { NaN, NaN, -0.89, 2.77, NaN, NaN, NaN, NaN, NaN, NaN },
					{ NaN, NaN, -0.27, 13.98, NaN, NaN, NaN, NaN, NaN, NaN }, { NaN, NaN, -0.83, 26.69, NaN, NaN, NaN, NaN, NaN, NaN },
					{ NaN, NaN, 1.84, 16.6, NaN, NaN, NaN, NaN, NaN, NaN },
					{ NaN, NaN, -0.6800000000000001, 20.79, NaN, NaN, NaN, NaN, NaN, NaN },
					{ NaN, NaN, NaN, NaN, 33.94, NaN, NaN, NaN, NaN, NaN }, { NaN, NaN, NaN, NaN, 20.85, NaN, NaN, NaN, NaN, NaN },
					{ NaN, NaN, NaN, NaN, 24.22, NaN, NaN, NaN, NaN, NaN }, { NaN, NaN, NaN, NaN, 22.97, NaN, NaN, NaN, NaN, NaN },
					{ NaN, NaN, NaN, NaN, 22.31, NaN, NaN, NaN, NaN, NaN }, { NaN, NaN, NaN, NaN, 14.87, NaN, NaN, NaN, NaN, NaN },
					{ NaN, NaN, NaN, NaN, 26.09, NaN, NaN, NaN, NaN, NaN }, { NaN, NaN, NaN, NaN, 23.99, NaN, NaN, NaN, NaN, NaN },
					{ NaN, NaN, NaN, NaN, 18.37, NaN, NaN, NaN, NaN, NaN }, { NaN, NaN, NaN, NaN, 21.41, 8.82, NaN, NaN, NaN, NaN },
					{ NaN, NaN, NaN, NaN, 30.19, 17.41, NaN, NaN, NaN, NaN }, { NaN, NaN, NaN, NaN, 19.45, 9.42, NaN, NaN, NaN, NaN },
					{ NaN, NaN, NaN, NaN, 14.93, 8.43, NaN, NaN, NaN, NaN }, { NaN, NaN, NaN, NaN, 16.69, 4.34, NaN, NaN, NaN, NaN },
					{ NaN, NaN, NaN, NaN, 8.359999999999999, 8.32, NaN, NaN, NaN, NaN },
					{ NaN, NaN, NaN, NaN, 26.02, 12.98, NaN, NaN, NaN, NaN }, { NaN, NaN, NaN, NaN, 10.13, 3.21, NaN, NaN, NaN, NaN },
					{ NaN, NaN, NaN, NaN, NaN, NaN, 12.92, NaN, NaN, NaN }, { NaN, NaN, NaN, NaN, NaN, NaN, 2.07, NaN, NaN, NaN },
					{ NaN, NaN, NaN, NaN, NaN, NaN, 10.77, NaN, NaN, NaN }, { NaN, NaN, NaN, NaN, NaN, NaN, 5.62, NaN, NaN, NaN },
					{ NaN, NaN, NaN, NaN, NaN, NaN, 3.85, NaN, NaN, NaN }, { NaN, NaN, NaN, NaN, NaN, NaN, 11.62, NaN, NaN, NaN },
					{ NaN, NaN, NaN, NaN, NaN, NaN, 12.58, NaN, NaN, NaN }, { NaN, NaN, NaN, NaN, NaN, NaN, 3.06, NaN, NaN, NaN },
					{ NaN, NaN, NaN, NaN, NaN, NaN, 18.84, NaN, NaN, NaN }, { NaN, NaN, NaN, NaN, NaN, NaN, NaN, -1.11, -0.7, 3.17 },
					{ NaN, NaN, NaN, NaN, NaN, NaN, NaN, 0.66, 0.13, 2.83 }, { NaN, NaN, NaN, NaN, NaN, NaN, NaN, 0.82, -0.17, 4.16 },
					{ NaN, NaN, NaN, NaN, NaN, NaN, NaN, -0.34, -0.23, -0.83 },
					{ NaN, NaN, NaN, NaN, NaN, NaN, NaN, 0.48, 0.5600000000000001, 2.24 },
					{ NaN, NaN, NaN, NaN, NaN, NaN, NaN, -0.25, 0.25, 7.56 }, { NaN, NaN, NaN, NaN, NaN, NaN, NaN, -0.71, -3.18, 0.73 },
					{ NaN, NaN, NaN, NaN, NaN, NaN, NaN, -1.08, -0.55, 5.58 }
			},
			{ { 536, 1236, 466, 1379, 919, 1111, 720, 458, 451, 519 }, { 369, 1032, 362, 1075, 711, 1168, 459, 336, 338, 424 },
					{ 792, 1687, 684, 1671, 1213, 1377, 897, 504, 448, 557 }, { 578, 1667, 472, 1646, 994, 1792, NaN, NaN, NaN, NaN },
					{ 653, 1802, 552, 1541, 1318, 1725, NaN, NaN, NaN, NaN }, { 769, 1185, 442, 1156, 901, 1676, 578, 365, 400, 434 },
					{ 1020, 1753, 487, 1475, 831, 1588, 612, 437, 422, 508 }, { 1366, 1607, 640, 1608, 1062, 1364, 843, 458, 576, 631 },
					{ 556, 1288, 526, 1353, 1251, 1726, NaN, NaN, NaN, NaN }, { NaN, NaN, NaN, NaN, NaN, NaN, 475, 339, 359, 449 },
					{ NaN, NaN, NaN, NaN, NaN, NaN, 1054, 615, 582, 676 }, { NaN, NaN, NaN, NaN, NaN, NaN, 988, 627, 478, 762 },
					{ NaN, NaN, NaN, NaN, NaN, NaN, 901, 479, 526, 751 }
			},
			{ { -0.0308, 0.1726, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN }, { 0.19, -0.0071, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN },
					{ -0.074, 0.067, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN }, { 0.0383, 0.2959, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN },
					{ -0.1171, 0.2379, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN }, { 0.0115, 0.0336, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN },
					{ 0.0494, 0.1544, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN }, { 0.0369, 0.3496, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN },
					{ 0.0187, 0.215, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN }, { 0.0524, 0.1846, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN },
					{ 0.0852, 0.1673, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN }, { 0.1068, 0.207, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN },
					{ 0.1217, 0.209, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN }, { 0.0123, 0.158, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN },
					{ -0.0475, 0.2306, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN },
					{ -0.0345, 0.1154, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN },
					{ -0.0808, 0.2935, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN }, { 0.0136, 0.24, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN },
					{ 0.041, 0.2966, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN }, { -0.0303, 0.1273, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN },
					{ 0.1712, 0.3195, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN }, { 0.1155, 0.1361, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN },
					{ 0.2328, 0.0814, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN }, { -0.0298, 0.2817, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN },
					{ -0.0447, 0.138, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN }, { 0.1699, 0.2192, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN },
					{ 0.1346, -0.0064, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN }, { 0.1755, 0.1453, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN },
					{ 0.1777, 0.2262, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN }, { -0.0164, 0.1258, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN },
					{ -0.0592, 0.06660000000000001, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN },
					{ 0.1593, 0.2792, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN }, { 0.0876, 0.1371, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN },
					{ 0.0315, 0.1489, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN }, { 0.014, 0.1134, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN },
					{ -0.0276, 0.2228, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN }, { 0.1831, 0.1894, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN },
					{ 0.0379, 0.2034, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN }, { 0.1281, 0.1765, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN },
					{ -0.0223, 0.2262, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN }, { 0.0633, 0.2747, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN },
					{ 0.0127, 0.1607, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN }, { 0.0945, 0.205, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN },
					{ 0.1182, 0.2808, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN }, { 0.1847, 0.2186, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN },
					{ 0.061, 0.218, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN }, { 0.0524, 0.2074, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN },
					{ 0.096, 0.2522, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN }, { 0.1241, 0.2214, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN },
					{ 0.1073, 0.1485, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN }, { 0.2181, 0.236, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN },
					{ 0.0078, 0.2305, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN }, { 0.0663, 0.1631, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN },
					{ 0.1585, 0.2308, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN }, { 0.0998, 0.0848, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN },
					{ 0.1765, 0.08350000000000001, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN },
					{ 0.0853, 0.2427, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN }, { 0.033, 0.1692, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN },
					{ 0.0815, 0.2203, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN }, { 0.1761, 0.1281, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN },
					{ 0.1023, 0.1428, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN },
					{ 0.09080000000000001, 0.2936, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN },
					{ 0.1116, 0.2347, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN }, { 0.2318, 0.1583, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN },
					{ 0.09810000000000001, 0.2383, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN },
					{ 0.0069, 0.088, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN }, { 0.1692, 0.2104, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN },
					{ 0.1098, 0.1932, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN }, { 0.0901, 0.2671, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN },
					{ 0.138, 0.1729, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN }, { 0.0854, 0.1271, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN },
					{ 0.0491, 0.1824, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN }, { NaN, 0.2966, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN },
					{ NaN, NaN, 0.0469, 0.125, NaN, NaN, NaN, NaN, NaN, NaN }, { NaN, NaN, 0.0417, 0.1823, NaN, NaN, NaN, NaN, NaN, NaN },
					{ NaN, NaN, 0.07290000000000001, 0.2604, NaN, NaN, NaN, NaN, NaN, NaN },
					{ NaN, NaN, 0.0625, 0.1615, NaN, NaN, NaN, NaN, NaN, NaN }, { NaN, NaN, 0.0625, -0.0052, NaN, NaN, NaN, NaN, NaN, NaN },
					{ NaN, NaN, -0.0469, 0.09900000000000001, NaN, NaN, NaN, NaN, NaN, NaN },
					{ NaN, NaN, -0.0104, 0.0573, NaN, NaN, NaN, NaN, NaN, NaN }, { NaN, NaN, 0.026, 0.1719, NaN, NaN, NaN, NaN, NaN, NaN },
					{ NaN, NaN, -0.2448, 0.3281, NaN, NaN, NaN, NaN, NaN, NaN }, { NaN, NaN, -0.026, 0.2708, NaN, NaN, NaN, NaN, NaN, NaN },
					{ NaN, NaN, -0.0469, 0.1927, NaN, NaN, NaN, NaN, NaN, NaN },
					{ NaN, NaN, -0.0625, 0.1302, NaN, NaN, NaN, NaN, NaN, NaN }, { NaN, NaN, 0.0573, 0.125, NaN, NaN, NaN, NaN, NaN, NaN },
					{ NaN, NaN, 0.0156, 0.224, NaN, NaN, NaN, NaN, NaN, NaN }, { NaN, NaN, -0.0469, 0.2135, NaN, NaN, NaN, NaN, NaN, NaN },
					{ NaN, NaN, 0.0417, 0.3594, NaN, NaN, NaN, NaN, NaN, NaN }, { NaN, NaN, 0.1406, 0.2135, NaN, NaN, NaN, NaN, NaN, NaN },
					{ NaN, NaN, 0.1666, 0.3021, NaN, NaN, NaN, NaN, NaN, NaN }, { NaN, NaN, 0.1823, 0.0365, NaN, NaN, NaN, NaN, NaN, NaN },
					{ NaN, NaN, 0.1093, 0.224, NaN, NaN, NaN, NaN, NaN, NaN }, { NaN, NaN, 0.0625, 0.3698, NaN, NaN, NaN, NaN, NaN, NaN },
					{ NaN, NaN, 0.2344, 0.1407, NaN, NaN, NaN, NaN, NaN, NaN }, { NaN, NaN, 0.2448, 0.3698, NaN, NaN, NaN, NaN, NaN, NaN },
					{ NaN, NaN, 0.1875, 0.2084, NaN, NaN, NaN, NaN, NaN, NaN },
					{ NaN, NaN, -0.09900000000000001, 0.276, NaN, NaN, NaN, NaN, NaN, NaN },
					{ NaN, NaN, -0.026, 0.1094, NaN, NaN, NaN, NaN, NaN, NaN }, { NaN, NaN, -0.1146, 0.0625, NaN, NaN, NaN, NaN, NaN, NaN },
					{ NaN, NaN, -0.0156, 0.0781, NaN, NaN, NaN, NaN, NaN, NaN },
					{ NaN, NaN, -0.0417, 0.2396, NaN, NaN, NaN, NaN, NaN, NaN }, { NaN, NaN, 0.2083, 0.1146, NaN, NaN, NaN, NaN, NaN, NaN },
					{ NaN, NaN, -0.1042, 0.026, NaN, NaN, NaN, NaN, NaN, NaN },
					{ NaN, NaN, 0.1562, 0.09900000000000001, NaN, NaN, NaN, NaN, NaN, NaN },
					{ NaN, NaN, 0.1667, 0.151, NaN, NaN, NaN, NaN, NaN, NaN }, { NaN, NaN, 0.1979, 0.2917, NaN, NaN, NaN, NaN, NaN, NaN },
					{ NaN, NaN, 0.0104, 0.1354, NaN, NaN, NaN, NaN, NaN, NaN }, { NaN, NaN, 0.1979, 0.1406, NaN, NaN, NaN, NaN, NaN, NaN },
					{ NaN, NaN, 0.1615, 0.1615, NaN, NaN, NaN, NaN, NaN, NaN }, { NaN, NaN, -0.0208, 0.2083, NaN, NaN, NaN, NaN, NaN, NaN },
					{ NaN, NaN, 0.0104, -0.1042, NaN, NaN, NaN, NaN, NaN, NaN },
					{ NaN, NaN, -0.0208, 0.1771, NaN, NaN, NaN, NaN, NaN, NaN }, { NaN, NaN, NaN, NaN, NaN, NaN, 0.1823, NaN, NaN, NaN },
					{ NaN, NaN, NaN, NaN, NaN, NaN, 0.09379999999999999, NaN, NaN, NaN },
					{ NaN, NaN, NaN, NaN, NaN, NaN, 0.1849, NaN, NaN, NaN }, { NaN, NaN, NaN, NaN, NaN, NaN, 0.1354, NaN, NaN, NaN },
					{ NaN, NaN, NaN, NaN, NaN, NaN, 0.2631, NaN, NaN, NaN }, { NaN, NaN, NaN, NaN, NaN, NaN, 0.2266, NaN, NaN, NaN },
					{ NaN, NaN, NaN, NaN, NaN, NaN, 0.1823, NaN, NaN, NaN }, { NaN, NaN, NaN, NaN, NaN, NaN, 0.3281, NaN, NaN, NaN },
					{ NaN, NaN, NaN, NaN, NaN, NaN, 0.4896, NaN, NaN, NaN }, { NaN, NaN, NaN, NaN, -0.1302, 0.0625, NaN, NaN, NaN, NaN },
					{ NaN, NaN, NaN, NaN, NaN, 0.1458, NaN, NaN, NaN, NaN }, { NaN, NaN, NaN, NaN, 0.0052, 0.1667, NaN, NaN, NaN, NaN },
					{ NaN, NaN, NaN, NaN, 0.0052, 0.1041, NaN, NaN, NaN, NaN },
					{ NaN, NaN, NaN, NaN, 0.09379999999999999, 0.09379999999999999, NaN, NaN, NaN, NaN },
					{ NaN, NaN, NaN, NaN, 0.0573, 0.2135, NaN, NaN, NaN, NaN },
					{ NaN, NaN, NaN, NaN, 0.0208, 0.09900000000000001, NaN, NaN, NaN, NaN },
					{ NaN, NaN, NaN, NaN, 0.0208, 0.1093, NaN, NaN, NaN, NaN }, { NaN, NaN, NaN, NaN, 0.1718, 0.1302, NaN, NaN, NaN, NaN },
					{ NaN, NaN, NaN, NaN, 0.1094, 0.224, NaN, NaN, NaN, NaN }, { NaN, NaN, NaN, NaN, 0.2344, 0.0469, NaN, NaN, NaN, NaN },
					{ NaN, NaN, NaN, NaN, 0.0677, 0.0469, NaN, NaN, NaN, NaN }, { NaN, NaN, NaN, NaN, NaN, 0.1198, NaN, NaN, NaN, NaN },
					{ NaN, NaN, NaN, NaN, 0.0521, 0.09379999999999999, NaN, NaN, NaN, NaN },
					{ NaN, NaN, NaN, NaN, 0.125, 0.2135, NaN, NaN, NaN, NaN }, { NaN, NaN, NaN, NaN, 0.2135, 0.1458, NaN, NaN, NaN, NaN },
					{ NaN, NaN, NaN, NaN, 0.1458, 0.1094, NaN, NaN, NaN, NaN }, { NaN, NaN, NaN, NaN, -0.0208, 0.0469, NaN, NaN, NaN, NaN },
					{ NaN, NaN, NaN, NaN, 0.1667, 0.0781, NaN, NaN, NaN, NaN },
					{ NaN, NaN, NaN, NaN, 0.0885, 0.07290000000000001, NaN, NaN, NaN, NaN },
					{ NaN, NaN, NaN, NaN, 0.1615, 0.1406, NaN, NaN, NaN, NaN }, { NaN, NaN, NaN, NaN, 0.2188, -0.0469, NaN, NaN, NaN, NaN },
					{ NaN, NaN, NaN, NaN, NaN, NaN, NaN, 0.0572, 0.2011, 0.2109 },
					{ NaN, NaN, NaN, NaN, NaN, NaN, NaN, 0.1459, 0.1154, 0.0848 },
					{ NaN, NaN, NaN, NaN, NaN, NaN, NaN, -0.0783, 0.2523, -0.0212 },
					{ NaN, NaN, NaN, NaN, NaN, NaN, NaN, 0.1348, 0.1221, -0.0276 },
					{ NaN, NaN, NaN, NaN, NaN, NaN, NaN, 0.0072, 0.1588, 0.119 },
					{ NaN, NaN, NaN, NaN, NaN, NaN, NaN, 0.2065, -0.046, 0.1567 },
					{ NaN, NaN, NaN, NaN, NaN, NaN, NaN, 0.1537, 0.0866, 0.18 },
					{ NaN, NaN, NaN, NaN, NaN, NaN, NaN, 0.0846, -0.021, 0.227 },
					{ NaN, NaN, NaN, NaN, NaN, NaN, NaN, 0.1826, 0.1245, -0.1694 },
					{ NaN, NaN, NaN, NaN, NaN, NaN, NaN, 0.0185, 0.0636, 0.1369 },
					{ NaN, NaN, NaN, NaN, NaN, NaN, NaN, 0.0999, 0.0228, 0.1998 },
					{ NaN, NaN, NaN, NaN, NaN, NaN, NaN, -0.0482, 0.2133, 0.2271 },
					{ NaN, NaN, NaN, NaN, NaN, NaN, NaN, 0.1724, 0.1535, 0.0385 },
					{ NaN, NaN, NaN, NaN, NaN, NaN, NaN, 0.025, 0.0776, 0.1151 },
					{ NaN, NaN, NaN, NaN, NaN, NaN, NaN, -0.0012, 0.1244, -0.0131 },
					{ NaN, NaN, NaN, NaN, NaN, NaN, NaN, 0.2506, 0.09039999999999999, 0.0019 },
					{ NaN, NaN, NaN, NaN, NaN, NaN, NaN, 0.1856, -0.0141, 0.2094 },
					{ NaN, NaN, NaN, NaN, NaN, NaN, NaN, 0.0179, 0.0279, 0.1796 },
					{ NaN, NaN, NaN, NaN, NaN, NaN, NaN, 0.0261, 0.0423, 0.3117 },
					{ NaN, NaN, NaN, NaN, NaN, NaN, NaN, 0.1333, -0.0239, 0.0202 },
					{ NaN, NaN, NaN, NaN, NaN, NaN, NaN, 0.1189, 0.1164, 0.0919 },
					{ NaN, NaN, NaN, NaN, NaN, NaN, NaN, 0.0861, 0.0161, 0.035 },
					{ NaN, NaN, NaN, NaN, NaN, NaN, NaN, 0.019, 0.07290000000000001, 0.2985 },
					{ NaN, NaN, NaN, NaN, NaN, NaN, NaN, 0.0185, 0.0258, 0.0394 },
					{ NaN, NaN, NaN, NaN, NaN, NaN, NaN, 0.0701, 0.0926, 0.2789 },
					{ NaN, NaN, NaN, NaN, NaN, NaN, NaN, -0.1294, -0.1472, 0.2721 },
					{ NaN, NaN, NaN, NaN, NaN, NaN, NaN, -0.0164, 0.0378, 0.1771 },
					{ NaN, NaN, NaN, NaN, NaN, NaN, NaN, 0.2169, 0.1895, 0.1817 },
					{ NaN, NaN, NaN, NaN, NaN, NaN, NaN, 0.3042, 0.238, -0.07729999999999999 },
					{ NaN, NaN, NaN, NaN, NaN, NaN, NaN, -0.1427, 0.2297, 0.1762 },
					{ NaN, NaN, NaN, NaN, NaN, NaN, NaN, 0.1957, 0.1307, 0.1276 },
					{ NaN, NaN, NaN, NaN, NaN, NaN, NaN, -0.0036, 0.0193, 0.2331 },
					{ NaN, NaN, NaN, NaN, NaN, NaN, NaN, 0.0375, 0.0954, 0.0451 },
					{ NaN, NaN, NaN, NaN, NaN, NaN, NaN, -0.0372, 0.1398, 0.1348 },
					{ NaN, NaN, NaN, NaN, NaN, NaN, NaN, -0.0235, -0.1278, 0.3947 },
					{ NaN, NaN, NaN, NaN, NaN, NaN, NaN, -0.0263, -0.0994, NaN }
			}
	};

	@Test
	public void statsTest() {
		GeneralModel gm = new GeneralModel();
		int var = 1;
		int s = 0;
		for (double[][] d : data) {
			for (double[] row : d) {
				gm.addData(s++, 1, var, row);
			}
			var++;
		}

		double[] expectedMeans = { 0.3275,
				20.0575,
				-0.1125,
				15.47125,
				20.87,
				9.11625,
				9.03666666666667,
				-0.19125,
				-0.48625,
				3.18 };

		CombinedStatsSTA[] stats = gm.calcStats();

		assertArrayEquals(expectedMeans, stats[0].getMeans().toArray(), 1e-13);

		CMRxGMProblemMaker maker = new CMRxGMProblemMaker();
		maker.setShrink(-1);
		maker.setModel(new double[][] { { 1 }, { 1 }, { 1 } });
		var = 1;
		for (double[][] d : data)
			maker.addCell(1, var++, d);

		CMRxSolver solver = new CMRxSolver();
		CMRSolution sol = solver.solve(maker.getProblem());

		//TODO: work out the correct value!
		//assertEquals(62.4968793544429, sol.getFStar(), 1e-8);
	}

	@Test
	public void covTest() {
		DoubleMatrix2D matrix = new DenseDoubleMatrix2D(data[0]);
		DoubleMatrix2D cov = MRUtil.covariance(matrix);

		DoubleMatrix2D expected = new DenseDoubleMatrix2D(new double[][] { { 2.64131875, -0.5867812499999991, 0, 0, 0, 0, 0, 0, 0, 0 },
				{ -0.5867812499999991, 33.86719375000001, 0, 0, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 1.70189375, 2.747528125000001, 0, 0, 0, 0, 0, 0 }, { 0, 0, 2.747528125000001, 47.81568593750001, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 41.65063529411766, 24.255003125, 0, 0, 0, 0 }, { 0, 0, 0, 0, 24.255003125, 17.8370234375, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 28.47211111111111, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0, 0, 0.5158109375, 0.3894671875, -0.05645 },
				{ 0, 0, 0, 0, 0, 0, 0, 0.3894671875, 1.1860234375, 0.9696 }, { 0, 0, 0, 0, 0, 0, 0, -0.05645, 0.9696, 6.1242 }
		});

		assertTrue(expected.equals(cov));

		ShrinkDiagonal sd = new ShrinkDiagonal(matrix);

		DoubleMatrix2D expectedReg = new DenseDoubleMatrix2D(new double[][] { { 2.64131875, -0.5867812499999991, 0, 0, 0, 0, 0, 0, 0, 0 },
				{ -0.5867812499999991, 33.86719375000001, 0, 0, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 1.70189375, 2.747528125000001, 0, 0, 0, 0, 0, 0 }, { 0, 0, 2.747528125000001, 47.81568593750001, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 41.65063529411766, 24.255003125, 0, 0, 0, 0 }, { 0, 0, 0, 0, 24.255003125, 17.8370234375, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 28.47211111111111, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0, 0, 0.5158109375, 0.3894671875, -0.05645 },
				{ 0, 0, 0, 0, 0, 0, 0, 0.3894671875, 1.1860234375, 0.9696 }, { 0, 0, 0, 0, 0, 0, 0, -0.05645, 0.9696, 6.1242 }
		});
		sd.shrink();

		assertTrue(expectedReg.equals(sd.getResult()));

		DoubleMatrix2D expectedW = new DenseDoubleMatrix2D(new double[][] {
				{ 3.040493048851036, 0.05267942555237896, 0, 0, 0, 0, 0, 0, 0, 0 },
				{ 0.05267942555237896, 0.2371295170912973, 0, 0, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 5.18128488818093, -0.2977208352196837, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, -0.2977208352196837, 0.1844163938104035, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 1.96114969029496, -1.254962536266125, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, -1.254962536266125, 2.155018766084008, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0, 0.3160987945318812, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 22.12617465372513, -8.537579131337132, 1.555641436423903 },
				{ 0, 0, 0, 0, 0, 0, 0, -8.537579131337132, 11.0423810100528, -1.826955189136732 },
				{ 0, 0, 0, 0, 0, 0, 0, 1.555641436423903, -1.826955189136732, 1.609880753481778 }
		});

		StatsSTA ssta = new StatsSTA(matrix);

		assertTrue(expectedW.equals(ssta.getWeights()));
		assertEquals(8.962696135046611, ssta.getLMValue(), 1e-14);
	}

	@Test
	public void fitsTest() {
		CMRxFitsGMProblemMaker pm = new CMRxFitsGMProblemMaker();
		pm.setShrink(-1);
		for (int i = 0; i < data.length; i++)
			pm.addCell(1, i, data[i]);
		pm.setModel(new double[] { 1, 1, 1 });

		//This is just a test to make sure we can do 1000 fast enough
		Fits sol = pm.solve(10, Runtime.getRuntime().availableProcessors());

		// System.out.println(sol.getP());
		// System.out.println(Arrays.toString(sol.getBadnesses()));
		// System.out.println(Arrays.toString(sol.getFits()));
		// System.out.println(Arrays.toString(sol.getTimes()));
	}
}
