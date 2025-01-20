public class VerletIntegration {

    // Constants
    static final double G = 4.3009e-3;       // gravitational constant in pc (km/s)^2 M_sun^-1
    static final double TIMESTEP = 3.156e15; // time step in seconds
    static final int SIZE = 5;               // number of masses
    static final int iteration = 10;

    // Arrays for masses, positions, and velocities
    static final double[] M     = {100, 200, 150, 300, 500}; // solar masses
    static final double[] X     = {10,  20,  30,  40,   50}; // pc (x positions)
    static final double[] Z     = {-30, 20,  10,  30,  -10}; // pc (initial z positions)
    static final double[] VZINIT= {0,   0,   0,   0,    0};  // km/s (initial velocities in z)

    // Background mass
    static double BM = 0.0;
    static {
        // sum(M) * 100
        double sumM = 0.0;
        for (double m : M) {
            sumM += m;
        }
        BM = sumM * 100;  // background mass
    }

    // Constants for background potential
    static final double a = 2700.0; // pc
    static final double b = 200.0;  // pc

    // Conversion factor (pc -> km), used in acceleration formula:
    //   1 pc = 3.0856e13 km
    static final double PC_TO_KM = 3.0856e13;


    // --------------------------------------------------------
    // Compute z-acceleration on the i-th mass from the j-th mass
    private static double accelFromMass(double xi, double xj, double zi, double zj, double mj) {
        /*
          Based on the gravitational potential: 
             phi_ij = - G * mj / r
          Acceleration on mass i in the z-direction is: 
             a_z = - d/dz_i (phi_ij)
                 = G * mj * (zj - zi) / r^3
          Then we divide by PC_TO_KM because positions are in pc,
          but we want acceleration in km/s^2.

          r = sqrt((xj - xi)^2 + (zj - zi)^2)
        */
        double dx = xj - xi;
        double dz = zj - zi;
        double r  = Math.sqrt(dx * dx + dz * dz);

        // Avoid division by zero if two masses overlap (not expected in normal usage)
        if (r < 1e-12) {
            return 0.0;
        }

        // a_z in km/s^2
        double az = G * mj * dz / (r * r * r);
        // Convert pc-based expression to correct units
        az = az / PC_TO_KM;

        return az;
    }


    // --------------------------------------------------------
    // Compute the z-acceleration due to the background mass
    // (i.e., partial derivative of the background potential).
    private static double accelFromBackground(double x, double z) {
        /*
           The background potential is:
               phi_b = - G * BM / sqrt( x^2 + ( sqrt(z^2 + b^2 ) + a )^2 )

           The z-acceleration is - d/dz(phi_b).  One can work it out
           explicitly, or replicate the symbolic result from your Python code.

           If we do it explicitly, the final form (in km/s^2) is typically:

               a_b_z = - G * BM * d/dz(1/R)      / PC_TO_KM
                      = - G * BM / (R^2) * dR/dz / PC_TO_KM

               where R = sqrt( x^2 + (sqrt(z^2 + b^2) + a)^2 ).
        */
        double w = Math.sqrt(z * z + b * b); // sqrt(z^2 + b^2)
        double R = Math.sqrt(x * x + (w + a) * (w + a));

        // to avoid tiny denominators:
        if (R < 1e-12) {
            return 0.0;
        }

        // dR/dz for R = sqrt( x^2 + (w + a)^2 ),  w = sqrt(z^2 + b^2 )
        //    dR/dz = ( (w + a)/R ) * (z / w )
        // Then the final acceleration a_b_z = - G * BM / R^2 * dR/dz, all / PC_TO_KM
        double dRdz = 0.0;
        if (w > 1e-12) {
            dRdz = ((w + a) / R) * (z / w);
        }
        double aBz = - G * BM / (R * R) * dRdz;

        // convert to km/s^2
        aBz = aBz / PC_TO_KM;

        return aBz;
    }


    // --------------------------------------------------------
    // Zaccel: total z-acceleration for each mass at given Z-list
    private static double[] Zaccel(double[] currentZ) {
        double[] acceleration = new double[SIZE];

        for (int i = 0; i < SIZE; i++) {
            // sum of pairwise accelerations from all other masses
            double sumAccel = 0.0;
            for (int j = 0; j < SIZE; j++) {
                if (j == i) continue;
                sumAccel += accelFromMass(X[i], X[j], currentZ[i], currentZ[j], M[j]);
            }
            // Add background acceleration
            sumAccel += accelFromBackground(X[i], currentZ[i]);

            acceleration[i] = sumAccel;
        }

        return acceleration;
    }


    // --------------------------------------------------------
    // Compute the initial positions after the first step
    // (from the usual Verlet "first step" formula).
    private static double[] firstStep() {
        double[] acceleration = Zaccel(Z); // acceleration at the initial Z
        double[] nextZ = new double[SIZE];

        // position(t+dt) = position(t) + v(t)*dt + 0.5*a(t)*dt^2
        // velocities are in km/s and distances in pc.
        // convert velocity from km/s to pc/s => v / PC_TO_KM.
        // multiply by TIMESTEP(s) to get pc. Similarly for a -> pc/s^2.
        for (int i = 0; i < SIZE; i++) {
            double zOld = Z[i];
            double vZ_km_s = VZINIT[i];
            double aZ_km_s2 = acceleration[i];

            // Convert velocity from km/s to pc/s
            double vZ_pc_s = vZ_km_s / PC_TO_KM;
            // Convert acceleration from km/s^2 to pc/s^2
            double aZ_pc_s2 = aZ_km_s2 / PC_TO_KM;

            // deltaZ = vZ_pc_s * TIMESTEP + 0.5 * aZ_pc_s2 * TIMESTEP^2
            double deltaZ = vZ_pc_s * TIMESTEP + 0.5 * aZ_pc_s2 * (TIMESTEP * TIMESTEP);

            nextZ[i] = zOld + deltaZ;
        }
        return nextZ;
    }


    // --------------------------------------------------------
    // Perform the full Verlet integration
    private static double[][] timeStep(double[] zPrev, double[] zCurr) {
        // We'll store all positions in an array of length (iteration+2).
        // each entry is an array [z1, z2, ..., zSIZE].
        double[][] positions = new double[iteration + 2][];
        
        // positions[0] -> initial
        positions[0] = zPrev.clone();
        // positions[1] -> after first step
        positions[1] = zCurr.clone();

        double[] currAccel = Zaccel(zCurr);

        for (int i = 0; i < iteration; i++) {
            double[] zNext = new double[SIZE];

            // zNext = 2*zCurr - zPrev + a(t)*dt^2 (in pc, so again convert acceleration)
            for (int j = 0; j < SIZE; j++) {
                double aZ_pc_s2 = currAccel[j] / PC_TO_KM; // convert to pc/s^2
                // 2zCurr - zPrev + a * dt^2
                zNext[j] = 2.0 * zCurr[j] - zPrev[j] + aZ_pc_s2 * (TIMESTEP * TIMESTEP);
            }

            // store
            positions[i + 2] = zNext.clone();
            
            // update for next iteration
            zPrev = zCurr;
            zCurr = zNext;
            currAccel = Zaccel(zNext);
        }

        return positions;
    }


    // --------------------------------------------------------
    // Main method to run
    public static void main(String[] args) {

        // Make time array similar to np.linspace(0, iteration*TIMESTEP, iteration+2)
        double[] timeTick = new double[iteration + 2];
        double maxTime = iteration * TIMESTEP;
        double step = maxTime / (iteration + 1);  // iteration+2 points
        for (int i = 0; i < iteration + 2; i++) {
            timeTick[i] = i * step;
        }

        // 1) Compute first-step positions
        double[] initZ = firstStep();

        // 2) Perform the Verlet integration
        double[][] positionArray = timeStep(Z, initZ);

        // Print out the results
        // positionArray[i] is the array of z-values (length = SIZE) at timeTick[i].
        // For demonstration:
        for (int i = 0; i < positionArray.length; i++) {
            System.out.printf("t = %.3e s:\t", timeTick[i]);
            for (int j = 0; j < SIZE; j++) {
                System.out.printf("Z%d=%.6f pc ", j, positionArray[i][j]);
            }
            System.out.println();
        }
    }
}
