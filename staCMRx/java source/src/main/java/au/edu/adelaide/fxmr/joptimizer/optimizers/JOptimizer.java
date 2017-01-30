/*
 * Copyright 2011-2014 JOptimizer
 *
 *   Licensed under the Apache License, Version 2.0 (the "License");
 *   you may not use this file except in compliance with the License.
 *   You may obtain a copy of the License at
 *
 *       http://www.apache.org/licenses/LICENSE-2.0
 *
 *   Unless required by applicable law or agreed to in writing, software
 *   distributed under the License is distributed on an "AS IS" BASIS,
 *   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *   See the License for the specific language governing permissions and
 *   limitations under the License.
 */
package au.edu.adelaide.fxmr.joptimizer.optimizers;



/**
 * Convex Optimizer.
 * 
 * The algorithm selection is implemented as a Chain of Responsibility pattern,
 * and this class is the client of the chain.
 * 
 * @see "S.Boyd and L.Vandenberghe, Convex Optimization"
 * @author <a href="mailto:alberto.trivellato@gmail.com">alberto trivellato</a>
 */
class JOptimizer {

	public static final int DEFAULT_MAX_ITERATION = 1000;
	public static final double DEFAULT_FEASIBILITY_TOLERANCE = 1.E-6;
	public static final double DEFAULT_TOLERANCE = 1.E-5;
	public static final double DEFAULT_TOLERANCE_INNER_STEP = 1.E-5;
	public static final double DEFAULT_KKT_TOLERANCE = 1.E-9;
	public static final double DEFAULT_ALPHA = 0.055;
	public static final double DEFAULT_BETA = 0.55;
	public static final double DEFAULT_MU = 10;
}
