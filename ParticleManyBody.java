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
     *argv[0] particle input file name
     *argv[1] parameter input file name
     *argv[2] output file name
     */
    public static void main (String[] argv) throws IOException {

	//Initiate scanners for the input and the output writer
	BufferedReader input = new BufferedReader(new FileReader(argv[0]));
	BufferedReader param = new BufferedReader(new FileReader(argv[1]));
	PrintWriter output = new PrintWriter(new FileWriter(argv[2]));
	Scanner parameters = new Scanner(param);
	Scanner in = new Scanner(input);

	//Scan the initial coordinates of the particles from the file
	Particle3D[] particleArray = new Particle3D[in.nextInt()];
	for (int i = 0; i < particleArray.length ; i++){
	    particleArray[i] = new Particle3D(in) ;
	}
	
	//Adjusting centre of mass
	//Total momentum and mass
	Vector3D totalMomentum = new Vector3D();
	double totalMass = 0.0;
	for (int i = 0; i < particleArray.length; i++){
	    totalMomentum.add(particleArray[i].getVelocity().
			 mult(particleArray[i].getMass()));
	    totalMass+=particleArray[i].getMass();
	}
       
	//Velocity of the centre of mass
	Vector3D comVelocity = totalMomentum.div(totalMass);
	System.out.printf("%s\n",comVelocity);
	//Correcting the velocities
	for (int i = 0; i < particleArray.length; i++){
	    particleArray[i].setVelocity(Vector3D.
					 subVector(particleArray[i].
						   getVelocity(),
						   comVelocity));
	}
       
	//Number of steps
	int numstep = parameters.nextInt(); 
	// Size of timestep
        double dt = parameters.nextDouble();     
        // Initial time
        double t = 0;
	//Gravitational Constant
	double g = parameters.nextDouble();

	//Double for total energy
	double e = totalEnergy(particleArray, g);
	//Double for maximum total energy
	double eMax = e;
	//Double for minimum total energy
	double eMin = e;
	
	//The start of the Verlet algorithm

	//Initiate the arrays for handling the forces
	Vector3D[] forceArray = new Vector3D[particleArray.length];
	for (int i = 0; i < forceArray.length; i++){
	    forceArray[i] = new Vector3D();
		}

	Vector3D[][] forceTable = new Vector3D[particleArray.length][particleArray.length];
	for (int i = 0; i < forceTable.length; i++){
	    for (int j = 0; j < forceTable.length; j++){
		forceTable[i][j] = new Vector3D();
		}
	}

	Vector3D[] newForceArray = new Vector3D[particleArray.length];
	for (int i = 0; i < newForceArray.length; i++){
	    newForceArray[i] = new Vector3D();
		}

	//Compute the initial force
	leapForceArray(particleArray, forceArray, forceTable, g);
 
        // Print the initial conditions to the files
	vmdEntry(particleArray, 1, output);
 
 
        //Loop over timesteps
        for (int i=0;i<numstep;i++){
 
             // Update the postion using current velocity
	    leapPositionArray(particleArray, forceArray, dt);
 
            // Update the force using current position
	    leapForceArray(particleArray, newForceArray, forceTable, g);
 
            // Update the velocity based on average of current and new force
            leapVelocityVerletArray(particleArray, forceArray,
				    newForceArray, dt);
	    
	    //Update the total energy using current position
	    e = totalEnergy(particleArray, g);
 	    //Check if current energy is minimum or maximum.
	    eMin = Math.min(eMin, e);
	    eMax = Math.max(eMax, e);
	    
	    //Update the old force
	    for(int j=0;j < particleArray.length; j++){
	     forceArray[j].copy(newForceArray[j]);
	    }
           
	    // Increase the time
            t = t + dt;
	    
	    // Print the current parameters to files
	    vmdEntry(particleArray, i+2, output);

	    /*    totalMomentum.setVector(0.0,0.0,0.0);
	totalMass = 0.0;
	for (int j = 0; j < particleArray.length; j++){
	    totalMomentum.add(particleArray[j].getVelocity().
			 mult(particleArray[j].getMass()));
	    totalMass+=particleArray[j].getMass();
	}
       
	//Velocity of the centre of mass
	comVelocity = totalMomentum.div(totalMass);
	System.out.printf("%s\n",comVelocity);
	//Correcting the velocities
	for (int j = 0; j < particleArray.length; j++){
	    particleArray[j].setVelocity(Vector3D.
					 subVector(particleArray[j].
						   getVelocity(),
						   comVelocity));
	}
	    */
        }

	//Print the maximum energy fluctuation
	System.out.printf("Maximum energy fluctuation: %10.7f\n", Math.abs(eMax-eMin));
	System.out.printf("Smaller than 1E-06? %s\n", 1.0E-06> Math.abs(eMax-eMin));
 
        // Close the output file
        output.close();
    }
    
    /**
     *Static method to calculate the gravitational force on a particle
     * 
     *@param  particle1 Particle3D for which the force has to be calculated
     *@param  particle2 Particle3D instance which acts on the first particle
     *@return Vector3D which is a gravitational force on particle1   
     */
    public static Vector3D getForce(Particle3D p1, Particle3D p2, double grav){
    double m1 = p1.getMass();
    double m2 = p2.getMass();
    Vector3D r = Particle3D.particleSeparation(p1,p2);
    return r.mult(-grav*m1*m2/(r.magSq()*r.mag()));
       }
 
    /**
     *Static method to calculate the total energy of a system E=KE+PE
     *
     *@param particle Particle3D for which the energy has to be calculated
     *@return double which is the potential energy of a particle 
     */
    public static double potentialEnergy(Particle3D p1, Particle3D p2, double grav){
	double m1 = p1.getMass();
	double m2 = p2.getMass();
	double r = Particle3D.particleSeparation(p1,p2).mag();
	return -grav*m1*m2/r;

    }
    public static double totalEnergy(Particle3D[] bodies, double grav){
	double energy = 0.0;
	for (int i = 0; i < bodies.length; i++){
	    energy+= bodies[i].kineticEnergy();
	    for (int j=i+1; j< bodies.length; j++){
		energy+=potentialEnergy(bodies[i], bodies[j], grav);
         }
	}
	return energy;
    }
    
    /**
     *Updates all positions of particles using the current velocities
     *
     *@param bodies an array of bodies to be integrated
     *@param oldforces old forces that acted on the bodies
     *@param forces forces that are currently acting on the bodies
     *@param dt time step for integration
     */
     public static void leapPositionArray(Particle3D[] bodies,
					  Vector3D[] forces, double dt){
	 for(int i=0;i<bodies.length;i++){
	     bodies[i].leapPosition(dt, forces[i]);
	 }
     }
      /**
     *Updates all velocities using the position of particles and forces by Velocity Verlet method
     *
     *@param bodies an array of bodies to be integrated
     *@param oldforces old forces that acted on the bodies
     *@param forces forces that are currently acting on the bodies
     *@param dt time step for integration
     */
     public static void leapVelocityVerletArray(Particle3D[] bodies,
						Vector3D[] oldForces,
						Vector3D[] forces,
						double dt){
	 for(int i=0;i<bodies.length;i++){
	     bodies[i].leapVelocity(dt,Vector3D.addVector(oldForces[i],
							  forces[i]).mult(0.5));

	 }
     }
      /**
     *Updates all forces using the position of particles
     *
     *@param bodies an array of bodies
     *@param oldforces old forces that acted on the bodies
     *@param forces forces that are currently acting on the bodies
     */
     public static void leapForceArray(Particle3D[] bodies, Vector3D[] forces, 
				       Vector3D[][] forcetable, double grav){
		 for(int i=0; i < bodies.length; i++){
		 for(int j=0; j < bodies.length; j++){
		     if (j>i){
			 forcetable[i][j] = getForce(bodies[i], bodies[j], grav); 
		     }
		     else if (j < i){
			 forcetable[i][j] = forcetable[j][i].mult(-1.0);
		     }
		 }
	 }
	 for(int i=0; i < bodies.length; i++){
	     forces[i].setVector(0.0, 0.0, 0.0);
	     for(int j=0; j < bodies.length; j++){
		 forces[i].add(forcetable[i][j]);
	      		     }
	     	 }
     }
      /**
     *Writes out particle’s parameters in format suitable for a VMD trajectory file
     *@param bodies a Particle3D array of bodies
     *@param step an integer representing the id of the step to write into the file
     *@param output Printwriter representing the output file 
     */
public static void vmdEntry(Particle3D[] bodies, int step, PrintWriter output){
	output.printf("%d\n",bodies.length);
	output.printf("Point = %d\n",step);
	for(int i=0; i<bodies.length;i++){
	    output.printf("%s\n", bodies[i]);
	}
    }    

}
   
