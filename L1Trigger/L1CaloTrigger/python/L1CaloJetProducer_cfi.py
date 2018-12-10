import FWCore.ParameterSet.Config as cms

L1CaloJetProducer = cms.EDProducer("L1CaloJetProducer",
    debug = cms.bool(False),
    HcalTpEtMin = cms.double(0.5),
    EcalTpEtMin = cms.double(0.5),
    EtMinForSeedHit = cms.double(2.5),
    l1CaloTowers = cms.InputTag("L1EGammaClusterEmuProducer","L1CaloTowerCollection","L1AlgoTest"),
    L1CrystalClustersInputTag = cms.InputTag("L1EGammaClusterEmuProducer", "L1EGXtalClusterEmulator", "L1AlgoTest"),

	emFractionBins = cms.vdouble([ 0.00,0.10,0.15,0.19,0.23,0.27,0.32,0.36,0.42,0.52,1.05]),
	absEtaBins = cms.vdouble([ 0.00,0.30,0.70,1.00,1.20,2.00]),
	jetPtBins = cms.vdouble([ 0.0,5.0,7.5,10.0,12.5,15.0,17.5,20.0,22.5,25.0,27.5,30.0,35.0,40.0,45.0,50.0,55.0,60.0,65.0,70.0,75.0,80.0,85.0,90.0,95.0,100.0,110.0,120.0,130.0,140.0,150.0,160.0,170.0,180.0,190.0,200.0,225.0,250.0,275.0,300.0,325.0,400.0,500.0]),
	jetCalibrations = cms.vdouble([
		4.789, 2.628, 2.038, 1.394, 1.220, 0.898, 0.775, 0.743, 0.628, 0.608, 0.545, 0.474, 0.438, 0.427, 0.421, 0.424, 0.412, 0.436, 0.440, 0.443, 0.438, 0.493, 0.516, 0.501, 0.485, 0.535, 0.563, 0.605, 0.575, 0.602, 0.608, 0.637, 0.642, 0.661, 0.692, 0.717, 0.775, 0.852, 0.839, 0.927, 1.002, 1.243,
		4.336, 2.391, 1.684, 1.339, 1.194, 0.960, 0.802, 0.709, 0.654, 0.573, 0.531, 0.468, 0.422, 0.389, 0.365, 0.344, 0.332, 0.344, 0.351, 0.347, 0.369, 0.369, 0.376, 0.405, 0.426, 0.449, 0.463, 0.470, 0.495, 0.482, 0.527, 0.524, 0.560, 0.571, 0.576, 0.594, 0.626, 0.663, 0.734, 0.789, 0.845, 1.016,
		4.512, 2.232, 1.709, 1.333, 1.083, 0.902, 0.768, 0.683, 0.588, 0.548, 0.499, 0.452, 0.393, 0.353, 0.325, 0.310, 0.293, 0.275, 0.282, 0.265, 0.277, 0.270, 0.274, 0.281, 0.316, 0.331, 0.352, 0.362, 0.394, 0.413, 0.417, 0.433, 0.443, 0.490, 0.481, 0.489, 0.515, 0.542, 0.584, 0.621, 0.661, 0.790,
		5.745, 2.897, 1.977, 1.653, 1.333, 0.965, 0.975, 0.815, 0.651, 0.613, 0.538, 0.496, 0.436, 0.382, 0.348, 0.319, 0.297, 0.285, 0.266, 0.262, 0.259, 0.268, 0.270, 0.270, 0.255, 0.290, 0.298, 0.324, 0.338, 0.355, 0.353, 0.369, 0.388, 0.366, 0.374, 0.417, 0.436, 0.467, 0.475, 0.489, 0.521, 0.619,
		5.968, 4.764, 3.933, 3.241, 2.809, 2.347, 2.312, 1.990, 1.820, 1.679, 1.529, 1.237, 1.183, 1.048, 0.995, 0.932, 0.877, 0.864, 0.829, 0.903, 0.858, 0.970, 0.943, 1.021, 1.109, 1.053, 1.180, 1.257, 1.366, 1.378, 1.473, 1.456, 1.429, 1.337, 1.320, 1.438, 1.213, 1.114, 1.068, 0.859, 0.603, 0.511,
		2.950, 2.626, 1.792, 1.375, 1.133, 0.938, 0.857, 0.746, 0.670, 0.623, 0.509, 0.512, 0.446, 0.422, 0.381, 0.385, 0.423, 0.408, 0.440, 0.472, 0.461, 0.513, 0.527, 0.480, 0.541, 0.553, 0.577, 0.606, 0.611, 0.616, 0.648, 0.647, 0.674, 0.700, 0.689, 0.725, 0.756, 0.783, 0.839, 0.842, 0.978, 1.078,
		2.910, 2.746, 1.725, 1.337, 1.031, 0.903, 0.760, 0.666, 0.637, 0.561, 0.519, 0.460, 0.417, 0.379, 0.334, 0.325, 0.335, 0.325, 0.350, 0.353, 0.365, 0.404, 0.395, 0.430, 0.460, 0.479, 0.496, 0.511, 0.526, 0.549, 0.569, 0.570, 0.592, 0.581, 0.592, 0.625, 0.647, 0.673, 0.697, 0.714, 0.825, 0.877,
		3.490, 2.472, 1.773, 1.321, 1.055, 0.892, 0.765, 0.669, 0.604, 0.510, 0.502, 0.432, 0.374, 0.330, 0.310, 0.298, 0.270, 0.272, 0.286, 0.298, 0.295, 0.319, 0.315, 0.344, 0.354, 0.389, 0.408, 0.421, 0.434, 0.453, 0.472, 0.479, 0.499, 0.484, 0.518, 0.522, 0.528, 0.550, 0.572, 0.592, 0.634, 0.713,
		5.233, 2.764, 1.881, 1.578, 1.370, 1.051, 0.935, 0.790, 0.688, 0.584, 0.533, 0.468, 0.419, 0.382, 0.310, 0.294, 0.274, 0.269, 0.254, 0.269, 0.260, 0.274, 0.291, 0.293, 0.316, 0.321, 0.348, 0.353, 0.370, 0.407, 0.386, 0.402, 0.424, 0.420, 0.434, 0.440, 0.457, 0.456, 0.490, 0.493, 0.521, 0.570,
		5.250, 4.461, 4.133, 3.243, 2.753, 2.246, 2.089, 1.756, 1.595, 1.410, 1.279, 1.136, 1.001, 0.853, 0.840, 0.757, 0.709, 0.817, 0.743, 0.735, 0.819, 0.837, 0.752, 0.717, 0.837, 0.841, 0.905, 0.888, 0.861, 1.012, 0.821, 0.885, 0.870, 0.872, 0.753, 0.747, 0.726, 0.634, 0.575, 0.496, 0.499, 0.455,
		4.170, 2.380, 1.712, 1.489, 1.108, 0.954, 0.806, 0.703, 0.653, 0.560, 0.505, 0.477, 0.428, 0.397, 0.382, 0.406, 0.429, 0.428, 0.437, 0.491, 0.505, 0.515, 0.576, 0.575, 0.569, 0.582, 0.592, 0.630, 0.631, 0.651, 0.675, 0.695, 0.676, 0.718, 0.716, 0.732, 0.780, 0.791, 0.842, 0.867, 0.916, 1.078,
		3.861, 2.313, 2.065, 1.357, 1.082, 0.903, 0.807, 0.673, 0.606, 0.513, 0.481, 0.444, 0.383, 0.354, 0.335, 0.333, 0.335, 0.354, 0.359, 0.410, 0.419, 0.460, 0.453, 0.486, 0.499, 0.512, 0.523, 0.542, 0.562, 0.572, 0.578, 0.602, 0.598, 0.619, 0.619, 0.650, 0.659, 0.683, 0.694, 0.715, 0.777, 0.860,
		4.564, 2.431, 1.746, 1.434, 1.034, 0.874, 0.773, 0.685, 0.566, 0.525, 0.475, 0.417, 0.370, 0.328, 0.299, 0.291, 0.291, 0.294, 0.301, 0.322, 0.351, 0.363, 0.393, 0.397, 0.423, 0.416, 0.446, 0.477, 0.469, 0.489, 0.486, 0.500, 0.511, 0.517, 0.527, 0.539, 0.538, 0.556, 0.571, 0.592, 0.629, 0.691,
		4.838, 2.826, 2.563, 1.653, 1.345, 1.001, 0.823, 0.730, 0.698, 0.587, 0.476, 0.455, 0.388, 0.328, 0.299, 0.277, 0.278, 0.290, 0.305, 0.313, 0.311, 0.313, 0.335, 0.354, 0.381, 0.379, 0.387, 0.425, 0.437, 0.422, 0.419, 0.424, 0.455, 0.455, 0.436, 0.466, 0.463, 0.489, 0.479, 0.504, 0.524, 0.569,
		5.115, 5.154, 3.434, 3.012, 2.607, 2.159, 2.103, 1.865, 1.419, 1.230, 1.129, 1.121, 0.915, 0.787, 0.697, 0.742, 0.682, 0.747, 0.748, 0.791, 0.809, 0.798, 0.822, 0.865, 0.898, 0.861, 0.748, 0.852, 0.814, 0.754, 0.675, 0.657, 0.677, 0.706, 0.596, 0.631, 0.611, 0.526, 0.507, 0.483, 0.487, 0.484,
		4.807, 2.462, 1.854, 1.442, 1.150, 0.918, 0.771, 0.708, 0.603, 0.525, 0.496, 0.439, 0.422, 0.385, 0.394, 0.432, 0.449, 0.492, 0.524, 0.525, 0.561, 0.583, 0.564, 0.587, 0.595, 0.620, 0.669, 0.657, 0.702, 0.677, 0.697, 0.685, 0.722, 0.710, 0.709, 0.738, 0.764, 0.776, 0.801, 0.818, 0.878, 0.979,
		4.667, 2.557, 1.866, 1.398, 1.121, 0.956, 0.739, 0.676, 0.584, 0.529, 0.479, 0.433, 0.378, 0.340, 0.343, 0.364, 0.382, 0.398, 0.439, 0.455, 0.471, 0.479, 0.493, 0.527, 0.528, 0.557, 0.558, 0.574, 0.590, 0.601, 0.620, 0.625, 0.622, 0.642, 0.625, 0.647, 0.675, 0.697, 0.704, 0.726, 0.760, 0.833,
		4.729, 2.440, 1.709, 1.351, 1.132, 0.901, 0.803, 0.667, 0.594, 0.547, 0.442, 0.401, 0.353, 0.313, 0.302, 0.314, 0.338, 0.354, 0.369, 0.397, 0.430, 0.411, 0.434, 0.450, 0.441, 0.473, 0.496, 0.481, 0.498, 0.522, 0.522, 0.503, 0.529, 0.521, 0.540, 0.555, 0.559, 0.570, 0.593, 0.611, 0.620, 0.689,
		5.426, 2.867, 2.101, 2.010, 1.395, 1.055, 0.877, 0.789, 0.637, 0.586, 0.567, 0.447, 0.360, 0.326, 0.317, 0.318, 0.321, 0.339, 0.357, 0.342, 0.385, 0.388, 0.378, 0.392, 0.382, 0.430, 0.428, 0.455, 0.447, 0.429, 0.453, 0.459, 0.458, 0.458, 0.459, 0.481, 0.484, 0.495, 0.514, 0.506, 0.535, 0.568,
		6.035, 4.869, 3.638, 3.066, 2.840, 2.098, 2.073, 1.827, 1.574, 1.547, 1.301, 1.123, 0.880, 0.721, 0.750, 0.735, 0.845, 0.890, 0.783, 0.881, 0.805, 0.816, 0.840, 0.958, 0.814, 0.790, 0.842, 0.770, 0.790, 0.836, 0.759, 0.751, 0.797, 0.665, 0.659, 0.648, 0.612, 0.552, 0.514, 0.536, 0.491, 0.478,
		4.194, 2.577, 1.948, 1.437, 1.122, 0.994, 0.841, 0.703, 0.614, 0.538, 0.505, 0.439, 0.409, 0.424, 0.444, 0.477, 0.518, 0.531, 0.562, 0.574, 0.596, 0.597, 0.646, 0.618, 0.617, 0.655, 0.677, 0.691, 0.688, 0.694, 0.710, 0.704, 0.717, 0.725, 0.762, 0.749, 0.747, 0.778, 0.809, 0.821, 0.860, 0.947,
		5.259, 3.061, 1.888, 1.402, 1.066, 0.908, 0.773, 0.662, 0.588, 0.530, 0.464, 0.431, 0.368, 0.357, 0.379, 0.418, 0.430, 0.462, 0.469, 0.494, 0.509, 0.554, 0.548, 0.549, 0.547, 0.585, 0.606, 0.598, 0.612, 0.610, 0.636, 0.629, 0.644, 0.647, 0.657, 0.671, 0.679, 0.691, 0.706, 0.728, 0.759, 0.837,
		5.977, 2.607, 1.945, 1.412, 1.096, 0.964, 0.767, 0.680, 0.586, 0.528, 0.447, 0.419, 0.354, 0.328, 0.360, 0.375, 0.406, 0.413, 0.428, 0.430, 0.466, 0.465, 0.495, 0.496, 0.484, 0.506, 0.509, 0.519, 0.512, 0.520, 0.551, 0.541, 0.555, 0.547, 0.536, 0.558, 0.574, 0.582, 0.597, 0.603, 0.629, 0.671,
		5.588, 3.444, 2.332, 1.594, 1.409, 0.937, 1.200, 0.744, 0.643, 0.578, 0.461, 0.411, 0.348, 0.378, 0.368, 0.369, 0.382, 0.406, 0.403, 0.436, 0.417, 0.456, 0.412, 0.492, 0.471, 0.436, 0.457, 0.476, 0.456, 0.479, 0.493, 0.483, 0.471, 0.477, 0.492, 0.498, 0.494, 0.502, 0.505, 0.526, 0.528, 0.571,
		6.428, 4.855, 4.248, 3.384, 2.735, 2.216, 1.855, 1.724, 1.756, 1.203, 1.150, 1.089, 0.921, 0.927, 0.915, 0.927, 0.854, 0.783, 0.876, 0.829, 0.999, 0.977, 0.787, 0.913, 0.737, 0.856, 0.856, 0.841, 0.731, 0.897, 0.717, 0.731, 0.807, 0.697, 0.724, 0.677, 0.562, 0.590, 0.535, 0.527, 0.497, 0.474,
		4.479, 2.721, 2.070, 1.478, 1.147, 0.937, 0.811, 0.711, 0.581, 0.531, 0.483, 0.433, 0.450, 0.462, 0.500, 0.529, 0.560, 0.591, 0.616, 0.629, 0.673, 0.662, 0.654, 0.638, 0.660, 0.703, 0.685, 0.713, 0.713, 0.718, 0.723, 0.757, 0.731, 0.753, 0.743, 0.768, 0.774, 0.790, 0.796, 0.833, 0.860, 0.929,
		4.374, 2.707, 1.856, 1.484, 1.141, 0.941, 0.755, 0.651, 0.583, 0.514, 0.466, 0.396, 0.397, 0.421, 0.443, 0.473, 0.494, 0.522, 0.531, 0.542, 0.574, 0.593, 0.586, 0.599, 0.594, 0.609, 0.615, 0.626, 0.617, 0.646, 0.638, 0.651, 0.657, 0.674, 0.661, 0.681, 0.685, 0.689, 0.722, 0.720, 0.755, 0.821,
		5.425, 2.612, 2.116, 1.482, 1.174, 0.966, 0.736, 0.650, 0.554, 0.520, 0.467, 0.387, 0.415, 0.414, 0.427, 0.432, 0.437, 0.462, 0.487, 0.475, 0.491, 0.494, 0.492, 0.504, 0.523, 0.533, 0.534, 0.522, 0.511, 0.546, 0.546, 0.558, 0.569, 0.558, 0.556, 0.574, 0.585, 0.576, 0.600, 0.609, 0.635, 0.665,
		5.135, 3.737, 2.561, 1.669, 1.385, 1.006, 0.925, 0.788, 0.678, 0.632, 0.511, 0.459, 0.426, 0.418, 0.440, 0.439, 0.438, 0.443, 0.454, 0.437, 0.428, 0.485, 0.436, 0.471, 0.481, 0.463, 0.485, 0.457, 0.484, 0.491, 0.516, 0.483, 0.484, 0.496, 0.481, 0.503, 0.502, 0.497, 0.520, 0.538, 0.533, 0.561,
		6.272, 4.733, 3.855, 3.186, 2.729, 2.206, 1.848, 1.873, 1.610, 1.324, 1.338, 1.201, 1.000, 0.929, 0.954, 0.982, 0.992, 0.877, 1.188, 1.042, 0.974, 0.821, 0.870, 0.876, 0.907, 0.865, 0.887, 0.909, 0.791, 0.795, 0.904, 0.739, 0.803, 0.717, 0.652, 0.682, 0.639, 0.612, 0.549, 0.532, 0.496, 0.489,
		5.076, 2.989, 2.162, 1.554, 1.156, 0.992, 0.784, 0.679, 0.589, 0.513, 0.479, 0.489, 0.493, 0.556, 0.575, 0.623, 0.599, 0.672, 0.682, 0.653, 0.704, 0.694, 0.718, 0.657, 0.721, 0.731, 0.750, 0.729, 0.755, 0.737, 0.750, 0.774, 0.764, 0.762, 0.781, 0.776, 0.782, 0.803, 0.803, 0.818, 0.847, 0.924,
		4.946, 3.019, 2.033, 1.536, 1.150, 0.923, 0.776, 0.677, 0.601, 0.498, 0.460, 0.444, 0.475, 0.498, 0.510, 0.528, 0.558, 0.569, 0.590, 0.596, 0.597, 0.596, 0.627, 0.628, 0.614, 0.628, 0.652, 0.653, 0.649, 0.660, 0.676, 0.666, 0.650, 0.660, 0.692, 0.688, 0.712, 0.713, 0.720, 0.729, 0.754, 0.800,
		4.427, 2.740, 2.053, 1.449, 1.225, 0.940, 0.820, 0.627, 0.603, 0.491, 0.449, 0.469, 0.462, 0.485, 0.496, 0.516, 0.502, 0.494, 0.505, 0.510, 0.523, 0.531, 0.558, 0.573, 0.564, 0.551, 0.540, 0.554, 0.555, 0.566, 0.555, 0.571, 0.566, 0.571, 0.597, 0.580, 0.581, 0.600, 0.606, 0.602, 0.625, 0.663,
		4.241, 3.483, 2.209, 1.491, 1.426, 0.986, 0.900, 0.699, 0.604, 0.535, 0.467, 0.491, 0.485, 0.474, 0.465, 0.506, 0.520, 0.487, 0.463, 0.478, 0.485, 0.468, 0.510, 0.472, 0.508, 0.490, 0.501, 0.473, 0.504, 0.496, 0.520, 0.480, 0.536, 0.506, 0.525, 0.503, 0.508, 0.516, 0.515, 0.529, 0.525, 0.556,
		7.473, 5.410, 4.072, 3.533, 3.058, 2.208, 1.952, 1.971, 1.456, 1.205, 1.239, 1.370, 1.260, 1.157, 0.987, 1.295, 1.084, 1.243, 1.189, 1.118, 1.034, 1.156, 0.885, 1.037, 0.996, 0.992, 0.922, 0.890, 0.876, 0.907, 0.817, 0.812, 0.770, 0.846, 0.808, 0.710, 0.709, 0.658, 0.657, 0.585, 0.538, 0.484,
		4.865, 3.154, 2.085, 1.572, 1.157, 0.973, 0.784, 0.679, 0.593, 0.536, 0.555, 0.576, 0.605, 0.625, 0.661, 0.686, 0.707, 0.718, 0.718, 0.709, 0.736, 0.735, 0.736, 0.717, 0.741, 0.741, 0.759, 0.768, 0.759, 0.772, 0.778, 0.784, 0.783, 0.759, 0.795, 0.774, 0.797, 0.796, 0.808, 0.827, 0.849, 0.900,
		4.737, 3.269, 2.239, 1.626, 1.249, 0.981, 0.790, 0.668, 0.562, 0.535, 0.534, 0.550, 0.569, 0.586, 0.593, 0.646, 0.613, 0.614, 0.647, 0.671, 0.654, 0.697, 0.658, 0.668, 0.670, 0.675, 0.669, 0.668, 0.676, 0.668, 0.688, 0.693, 0.703, 0.686, 0.708, 0.702, 0.720, 0.716, 0.714, 0.726, 0.744, 0.787,
		4.889, 3.227, 2.274, 1.559, 1.289, 0.984, 0.850, 0.709, 0.602, 0.598, 0.553, 0.553, 0.572, 0.570, 0.557, 0.571, 0.603, 0.533, 0.591, 0.557, 0.587, 0.611, 0.565, 0.603, 0.612, 0.597, 0.614, 0.579, 0.588, 0.588, 0.600, 0.580, 0.578, 0.592, 0.599, 0.604, 0.587, 0.594, 0.587, 0.621, 0.620, 0.656,
		5.137, 4.193, 3.042, 2.535, 1.455, 1.257, 0.869, 0.790, 0.639, 0.683, 0.670, 0.616, 0.493, 0.650, 0.552, 0.531, 0.561, 0.618, 0.557, 0.534, 0.559, 0.549, 0.518, 0.496, 0.571, 0.552, 0.500, 0.525, 0.515, 0.493, 0.495, 0.514, 0.488, 0.544, 0.512, 0.494, 0.506, 0.499, 0.523, 0.524, 0.525, 0.540,
		6.625, 5.453, 4.619, 3.430, 2.838, 2.788, 2.071, 1.695, 1.585, 1.590, 1.458, 1.563, 1.470, 1.242, 1.415, 1.362, 1.189, 1.240, 1.375, 1.296, 1.229, 1.032, 1.052, 1.132, 1.007, 1.104, 1.105, 0.997, 0.988, 0.905, 0.997, 0.819, 0.873, 0.840, 0.893, 0.773, 0.688, 0.640, 0.634, 0.677, 0.577, 0.552,
		4.997, 3.774, 2.310, 1.722, 1.236, 1.004, 0.818, 0.686, 0.678, 0.730, 0.689, 0.730, 0.747, 0.763, 0.800, 0.788, 0.814, 0.774, 0.803, 0.792, 0.790, 0.787, 0.808, 0.809, 0.820, 0.788, 0.795, 0.800, 0.783, 0.795, 0.806, 0.805, 0.818, 0.792, 0.779, 0.825, 0.815, 0.804, 0.823, 0.828, 0.835, 0.882,
		5.103, 3.753, 2.324, 1.687, 1.298, 1.010, 0.889, 0.691, 0.635, 0.729, 0.697, 0.710, 0.758, 0.687, 0.724, 0.758, 0.706, 0.722, 0.714, 0.711, 0.726, 0.707, 0.705, 0.719, 0.692, 0.719, 0.718, 0.695, 0.708, 0.730, 0.726, 0.703, 0.720, 0.714, 0.694, 0.712, 0.719, 0.709, 0.722, 0.723, 0.738, 0.775,
		5.070, 3.581, 2.159, 1.872, 1.482, 1.014, 0.799, 0.800, 0.745, 0.729, 0.792, 0.759, 0.780, 0.669, 0.669, 0.745, 0.638, 0.659, 0.642, 0.619, 0.618, 0.621, 0.599, 0.639, 0.598, 0.644, 0.641, 0.604, 0.621, 0.607, 0.605, 0.602, 0.575, 0.586, 0.589, 0.599, 0.586, 0.596, 0.596, 0.623, 0.610, 0.638,
		5.183, 4.639, 2.864, 1.899, 1.374, 1.329, 1.031, 0.848, 0.743, 0.737, 0.758, 0.842, 0.861, 0.733, 0.620, 0.683, 0.664, 0.626, 0.594, 0.637, 0.602, 0.549, 0.608, 0.600, 0.565, 0.542, 0.532, 0.533, 0.556, 0.538, 0.517, 0.472, 0.527, 0.521, 0.513, 0.524, 0.520, 0.512, 0.497, 0.510, 0.512, 0.557,
		7.907, 5.746, 4.750, 3.741, 3.428, 2.712, 2.096, 2.058, 2.515, 1.952, 2.099, 1.777, 1.779, 1.687, 1.550, 1.879, 1.503, 1.338, 1.240, 1.289, 0.982, 1.407, 1.230, 1.452, 1.370, 1.281, 1.254, 1.161, 1.115, 1.013, 1.084, 1.166, 0.934, 1.017, 1.005, 0.989, 0.906, 0.855, 0.810, 0.652, 0.680, 0.610,
		5.198, 4.274, 3.276, 2.340, 1.805, 1.585, 1.382, 1.455, 1.444, 1.419, 1.377, 1.331, 1.307, 1.321, 1.272, 1.170, 1.122, 1.092, 1.073, 1.131, 1.063, 1.064, 1.000, 1.010, 0.935, 0.965, 0.909, 0.917, 0.840, 0.850, 0.874, 0.867, 0.869, 0.868, 0.802, 0.822, 0.844, 0.846, 0.822, 0.893, 0.850, 0.915,
		2.600, 4.344, 3.310, 2.420, 1.932, 1.605, 1.557, 1.493, 1.537, 1.353, 1.293, 1.314, 1.224, 1.183, 1.158, 1.097, 1.024, 1.000, 0.971, 0.987, 0.940, 0.890, 0.858, 0.899, 0.865, 0.830, 0.819, 0.773, 0.737, 0.758, 0.742, 0.705, 0.733, 0.734, 0.724, 0.715, 0.744, 0.751, 0.680, 0.708, 0.779, 0.787,
		4.050, 4.329, 3.270, 2.536, 2.167, 1.808, 1.674, 1.661, 1.627, 1.432, 1.337, 1.417, 1.437, 1.288, 1.119, 0.971, 1.013, 0.982, 0.940, 0.866, 0.778, 0.862, 0.755, 0.785, 0.692, 0.732, 0.690, 0.715, 0.730, 0.622, 0.697, 0.615, 0.582, 0.601, 0.624, 0.675, 0.654, 0.651, 0.668, 0.667, 0.652, 0.707,
		5.450, 4.012, 3.714, 2.572, 2.316, 1.815, 1.730, 1.571, 1.489, 1.125, 1.330, 1.791, 1.297, 1.190, 1.107, 1.029, 0.951, 0.814, 0.812, 0.813, 0.826, 0.675, 0.804, 0.743, 0.708, 0.662, 0.550, 0.576, 0.543, 0.625, 0.630, 0.633, 0.717, 0.722, 0.506, 0.651, 0.537, 0.606, 0.671, 0.582, 0.676, 0.652,
		0.850, 4.931, 4.880, 4.171, 4.011, 3.674, 3.448, 4.030, 3.932, 3.782, 3.022, 3.018, 3.185, 3.331, 2.636, 2.730, 2.782, 2.522, 2.182, 2.581, 2.503, 2.125, 2.295, 2.289, 1.901, 2.025, 1.815, 1.913, 2.154, 1.891, 1.762, 1.660, 1.632, 1.409, 1.547, 1.479, 1.459, 1.313, 1.297, 1.374, 1.215, 1.077
	])
)

