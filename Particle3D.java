/**
 * Computer Modelling, Exercise 3: Particle3D class to simulate a particle with its
 * position, velocity, mass and label. Complete with getters, setters, constructors
 * and methods to time integration of a particle motion.
 *
 * @author Rapolas Daugintis
 * @author Rokas Zemaitis
 */
import java.io.*;
import java.util.Scanner;
public class Particle3D {

    /* ******************************************
     * Properties
     ********************************************/
    private double mass;
    private Vector3D position;
    private Vector3D velocity;
    private String label;

    // Getters
    
    /** Get the position of a particle.
     *
     * @return a Vector3D representing the position.
     */
    public Vector3D getPosition() { return new Vector3D(position); }

    /** Get the velocity of a particle.
     *
     * @return a Vector3D representing the velocity.
     */
    public Vector3D getVelocity() { return new Vector3D(velocity); }

    /** Get the mass of a particle.
     *
     * @return a double representing the mass.
     */
    public double getMass()     { return mass; }

    /**Get the label of a particle.
     *
     *@return a String representing the label of a particle.
     */
    public String getLabel() {return label;}
   
    //Setters

    /** Set the position of a particle.
     *
     * @param p a Vector3D representing the position.
     */
    public void setPosition(Vector3D p) {position = new Vector3D(p); }
    
    /** Set the velocity of a particle.
     *
     * @param v a Vector3D representing the velocity.
     */
    public void setVelocity(Vector3D v) {velocity = new Vector3D(v); }
    
    /** Set the mass of a particle.
     *
     * @param m a double representing the mass.
     */
    public void setMass(double m) {mass = m; }

    /** Set the label of a particle
     *
     *@param l a string representing the label of a particle.
     */
    public void setLabel(String l) {label = l; }

     /* ******************************************
     * Constructors
     ********************************************/
    
    /** Default constructor. Sets all properties to zero.
     */
    public Particle3D() {
	mass = 0.0;
        position = new Vector3D();
        velocity = new Vector3D();
	label = "Untitled";
    }

    /** Explicit constructor. Constructs a new Particle3D with
     * explicitly given position, velocity, mass and label.
     *
     * @param m a double that defines the mass.
     * @param p a Vector3D that defines the position.
     * @param v a Vector3D that defines the velocity.
     * @param l a String that is the label of the particle.
     */
    public Particle3D(double m, Vector3D p, Vector3D v, String l) {
        mass = m;
        position = new Vector3D(p);
        velocity = new Vector3D(v);
	label = l;
    }
    
     /**Constructs the Particle3D with values scannes from the input file
     *
     *@param scan Scanner which contains the input file
     */

    public Particle3D(Scanner scan) throws IOException{
	label = scan.next();
	mass = scan.nextDouble();
	position = new Vector3D(scan.nextDouble(),scan.nextDouble(),scan.nextDouble());
	velocity = new Vector3D(scan.nextDouble(),scan.nextDouble(),scan.nextDouble());
	
    }

    /* ******************************************
     * toString Method
     ********************************************/
    
    /** Returns a String representation of Particle3D.
     * Used to print a Particle3D instance.
     *
     *@return String representing the label and position of the particle.
     */
    public String toString() {
        return String.format("%s %.5f %.5f %.5f",getLabel(),position.getX(),
			     position.getY(),position.getZ());
    }

     /* ******************************************
     * Instance and Static Methods
     ********************************************/
     /**Returns kinetic energy of a particle 0.5*m*|v|^2.
     *
     *@return a double, which represents kinetic energy.
     */
    public double kineticEnergy(){
	return 0.5*mass*velocity.magSq();
	    }

   /** Time integration support: evolve the velocity
     * according to v(t+dt) = v(t) + f/m * dt.
     *
     * @param dt a double that is the timestep.
     * @param f a Vector3D that is the current force on the particle.
     */
    public void leapVelocity(double dt, Vector3D f) {
        velocity = Vector3D.addVector(velocity, f.mult(dt/mass));
    }
    
    /** Time integration support: evolve the position
     * according to r(t+dt) = r(t)+ v(t) * dt.
     *
     * @param dt a double that is the timestep.
     */
    public void leapPosition(double dt) {
        position = Vector3D.addVector(position, velocity.mult(dt));
    }

    /** Time integration support: evolve the position
     * according to r(t+dt) = r(t)+v(t) * dt + 0.5 * f(t)/m * dt^2.
     *
     * @param dt a double that is the timestep.
     * @param force a Vector3D that is the current force.
     */
    public void leapPosition(double dt, Vector3D force) {
        position = Vector3D.addVector(position,
				      Vector3D.addVector(velocity.mult(dt),
							 force.mult(dt*dt/(2.0*mass))));
    }
    /**Static method, which returns the separation between two Particles.
     *
     *@param a First Particle
     *@param b Second Particle
     *@return Vector3D representing the separation between a and b (a-b)
     */
    public static Vector3D particleSeparation(Particle3D a, Particle3D b){
	return Vector3D.subVector(a.getPosition(), b.getPosition());
    }
}
