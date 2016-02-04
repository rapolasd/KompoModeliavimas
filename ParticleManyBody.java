/**
 * The code to simulate an N-body Solar system using Verlet time
 * integration.
 *
 * @author Rapolas Daugintis
 * @author Rokas Zemaitis
 */
 
import java.io.*;
import java.util.Scanner;
public class ParticleManyBody {
 
    /**
     *Main method of the class, which scans the particle from the
     *input file, defined by first argument of the class and outputs
     *modelled trajectory.
     *
     *argv[0] input file name
     *argv[1] size of a timestep
     */
 
    public static void main (String[] argv) throws IOException {
 
    //Create new particle
    Particle3D p = new Particle3D();
     
    //Open the file which the values will be readed from
    BufferedReader file = new BufferedReader(new FileReader(argv[0]));
     
    //Attach a Scanner to the file to parse the numbers
    Scanner scan = new Scanner(file);
    //Scan the values from the file and assign them to Particle3D
    p.scan(scan);
 
    //Print the particle to the terminal
    System.out.printf("Particle: %s\n", p);
 
    // Open the output files for different plots
        PrintWriter outputX = new PrintWriter(new FileWriter("XvsT_Verlet"));
    PrintWriter outputY = new PrintWriter(new FileWriter("YvsT_Verlet"));
    PrintWriter outputE = new PrintWriter(new FileWriter("EvsT_Verlet"));
    PrintWriter outputXY = new PrintWriter(new FileWriter("YvsX_Verlet"));
 
    //Initial Force
    Vector3D force = getForce(p);
 
    //Initial  energy
    double E = getEnergy(p);
 
    /*Creating 2 variables for maximum and minimum energy to see
     *maximum energy fluctuation.
     */
    double Emin = E;
    double Emax = E;
 
    // Size of timestep
        double dt = Double.parseDouble(argv[1]);
    // Number of timesteps for 3 full rotations
        int numstep = (int)(6.0*Math.PI/dt);     
        // Initial time
        double t = 0;
 
    //The start of the Verlet algorithm
 
        // Print the initial conditions to the files
    outputX.printf("%f %f\n", t, p.getPosition().getX());
    outputY.printf("%f %f\n", t, p.getPosition().getY());
    outputE.printf("%f %f\n", t, E);
    outputXY.printf("%f %f\n",p.getPosition().getX(), p.getPosition().getY());
 
 
        //Loop over timesteps
        for (int i=0;i<numstep;i++){
 
             // Update the postion using current velocity
        p.leapPosition(dt, force);
 
            // Update the force using current position
        Vector3D forceNew = new Vector3D(getForce(p));
 
            // Update the velocity based on average of current and new force
            p.leapVelocity(dt,Vector3D.addVector(force, forceNew).mult(0.5));
 
        //Update the total energy using current position
        E = getEnergy(p);
 
        //Check if current energy is minimum or maximum.
        Emin = Math.min(Emin, E);
        Emax = Math.max(Emax, E);
 
        //Update the old force
        force = new Vector3D(forceNew);
 
            // Increase the time
            t = t + dt;
 
            // Print the current parameters to files
        outputX.printf("%f %f\n", t, p.getPosition().getX());
        outputY.printf("%f %f\n", t, p.getPosition().getY());
        outputE.printf("%f %f\n", t, E);
        outputXY.printf("%f %f\n",p.getPosition().getX(), p.getPosition().getY());
        }
    //Print the maximum energy fluctuation
    System.out.printf("Maximum energy fluctuation: %10.7f\n", Math.abs(Emax-Emin));
    System.out.printf("Smaller than 1E-06? %s\n", 1.0E-06> Math.abs(Emax-Emin));
 
        // Close the output files
        outputX.close();
    outputY.close();
    outputE.close();
    outputXY.close();
    }
    
    /**
     *Static method to calculate the gravitational force on a particle
     * 
     *@param  particle Particle3D for which the force has to be calculated
     *@return Vector3D which is a gravitational force on a particle   
     */
    public static Vector3D getForce(Particle3D particle){
    double m1 = particle.getMass();
    //Hard-wired second mass
    double m2 = 1;
    Vector3D r = new Vector3D(particle.getPosition());
        return r.mult(-m1*m2/(r.magSq()*r.mag()));
    }
 
    /**
     *Static method to calculate the total energy of a system E=KE+PE
     *
     *@param particle Particle3D for which the energy has to be calculated
     *@return double which is the potential energy of a particle 
     */
    public static double potentialEnergy(Particle3D a, Particle3D b){

    }
    
    /**
     *Updates all positions of particles using the current velocities
     *
     *@param bodies an array of bodies to be integrated
     *@param oldforces old forces that acted on the bodies
     *@param forces forces that are currently acting on the bodies
     *@param dt time step for integration
     */
     public static void leapPositionArray(Particle3D[] bodies, Vector3D[] oldforces,
     Vector3D[] forces, double dt){
     }
      /**
     *Updates all velocities using the position of particles and forces by Velocity Verlet method
     *
     *@param bodies an array of bodies to be integrated
     *@param oldforces old forces that acted on the bodies
     *@param forces forces that are currently acting on the bodies
     *@param dt time step for integration
     */
     public static void leapVelocityVerletArray(Particle3D[] bodies, Vector3D[] oldforces,
     Vector3D[] forces, double dt){
      
      
     }
      /**
     *Updates all forces using the position of particles
     *
     *@param bodies an array of bodies
     *@param oldforces old forces that acted on the bodies
     *@param forces forces that are currently acting on the bodies
     */
     public static void leapForceArray(Particle3D[] bodies, Vector3D[] oldforces,
     Vector3D[] forces, double dt){
     
     
     }
      /**
     *Writes out particle’s parameters in format suitable for a VMD trajectory file
     *@param b an array of bodies
     *@return Vector3D which is the total energy of a particle 
     */
      public static String vmdEntry(Particle3D[] b){
       
      }
      
     }
     
}
