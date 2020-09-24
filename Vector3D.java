/**
 * Computer Modelling, Exercise 2:
 * A class for 3D Vectors, complete with constructors, setters, getters,
 * static methods to calculate vector and scalar products,
 * instance methods for vector coordinates manipulation,
 * and a boolean operation for vector equality checking.
 *
 * @author Rokas Zemaitis s1307874
 * @author Rapolas Daugintis s1348455
 * @version "10/2015"
 *
 */
 
public class Vector3D {
 
    /*
     * Properties
     *
     */
    private double x;
    private double y;
    private double z;
 
    /*
     * Constructors
     */
 
    /**
     * Default constructor. Constructs a new Vector3D, with uninitialised
     * vector elements.
     */
 
    public Vector3D(){
    //Set as a zero vector
    this.setVector(0.0, 0.0, 0.0);
    }
 
    /**
     * Copy constructor. Constructs a new Vector3D by copying the coordinates
     * of another Vector3D instance.
     *
     * @param original the Vector3D to be copied
     */
 
    public Vector3D (Vector3D original) {
    setVector(original.getX(), original.getY(), original.getZ());
        }
    /**
     * Explicit constructor. Constructs a new Vector3D from explicitly given
     * x, y and z elements.
     * 
     * @param xx a double giving the x coordinate of the new Vector3D
     * @param yy a double giving the y coordinate of the new Vector3D
     * @param zz a double giving the z coordinate of the new Vector3D
     */
 
    public Vector3D(double xx, double yy, double zz){
    setX(xx);
    setY(yy);
    setZ(zz);
    }
 
    /*
     * Setters and getters
     */
 
    /**
     * Set method to set all the coordinates at once.
     *
     * @param xx a double to set the x coordinate
     * @param yy a double to set the y coordinate
     * @param zz a double to sat the z coordinate
     */
    public void setVector(double xx, double yy, double zz) {
    x=xx;
    y=yy;
    z=zz;
    }
 
    //Setters
 
    /** Sets the x element only.
     *
     * @param xx a double to set the x element
     */
    public void setX(double xx){ x = xx; }
 
   /** Sets the y element only.
     *
     * @param yy a double to set the y element
     */
    public void setY(double yy){ y = yy; }
 
   /** Sets the z element only.
     *
     * @param zz a double to set the z element
     */
    public void setZ(double zz){ z = zz; }
 
    //Getters
 
    /**
     * Gets the x element of Vector3D
     *
     * @return a double instance representing the Vector3D' x element.
     */
    public double getX(){ return x;}
 
    /**
     * Gets the y element of Vector3D
     *
     * @return a double instance representing the Vector3D' y element.
     */
    public double getY(){ return y;}
 
    /**
     * Gets the z element of Vector3D
     *
     * @return a double instance representing the Vector3D' z element.
     */
    public double getZ(){ return z;}
 
    /**
     * Returns a String representation of Vector3D.<br>
     * For example, if a Vector3D has parameters x as 1.0,
     * y as- 2.0 and z as 3.0, the output is
     * "(1.0, -2.0, 3.0)".
     *
     * @return a string representation of the Vector3D instance
     */
    public String toString() {
        return  "(" + getX() +", " + getY() + ", " + getZ() + ")";
    }
 
    /*
     * Instance methods
     *
     *
     */
 
    /** Calculates the magnitude of the Vector3D squared.
     *
     * @return a double representing the magnitude of the Vector3D squared.
     */
 
    public double magSq(){
    return this.getX()*this.getX()+
        this.getY()*this.getY()+
        this.getZ()*this.getZ();
    }
 
 
    /** Calculates the magnitude of the Vector3D.
     *
     * @return a double representing the magnitude of the Vector3D.
     */
 
    public double mag(){
    return Math.sqrt(this.magSq());
    }
 

    /** Multiplies a Vector3D by a double.
     *
     * @param a the double being multiplied by a vector
     * @return a Vector3D multiplied by a double.
     */
 
    public Vector3D mult(double a){
    return new Vector3D(this.getX()*a,
                this.getY()*a,
                this.getZ()*a);
    }
 
    /** Divides a Vector3D by a double.
     *
     * @param a the double being divided by a vector
     * @return a Vector3D divided by a double.
     */
 
    public Vector3D div(double a){
    return new Vector3D(this.getX()/a,
                this.getY()/a,
                this.getZ()/a);
    }
     /** Adds a vector, changing the values of the vector that invokes this method.
     *
     * @param a the vector that is added to
     */
    public void add(Vector3D a){
    x+=a.getX();
    y+=a.getY();
    z+=a.getZ();
    }
    /** Subtracts a vector from the vector that invokes this method.
      *
     * @param a the vector that subtracts.
     */
    public void sub(Vector3D a){
    x-=a.getX();
    y-=a.getY();
    z-=a.getZ();
    }

    /** Copies coordinates from another Vector3D instance
     *
     * @param original the Vector3D instance to be copied from
     */
    public void copy(Vector3D original){
	this.setX(original.getX());
	this.setY(original.getY());
	this.setZ(original.getZ());
	
    }
    /*
     * Static methods
     *
     *
     */
 
    /** Adds two vectors together.
     *
     * @param a the first Vector3D
     * @param b the second Vector3D
     * @return a Vector3D instance representing the sum of the two vectors.
     */
 
    public static Vector3D addVector(Vector3D a, Vector3D b){
    return new Vector3D(a.getX()+b.getX(),
                a.getY()+b.getY(),
                a.getZ()+b.getZ());
    }
     
    /** Subtracts a vector from another vector.
     *
     * @param a the first Vector3D
     * @param b the second Vector3D
     * @return a Vector3D instance representing first Vector3D subtracted by a second Vector3D.
     */
 
    public static Vector3D subVector(Vector3D a, Vector3D b){
    return new Vector3D(a.getX()-b.getX(),
                a.getY()-b.getY(),
                a.getZ()-b.getZ());
    }
 
    /** Calculates the dot product of two vectors.
     *
     * @param a the first Vector3D
     * @param b the second Vector3D
     * @return a double representing the dot product of the two vectors.
     */
 
    public static double dotVector(Vector3D a, Vector3D b){
    return (a.getX()*b.getX()+
        a.getY()*b.getY()+
        a.getZ()*b.getZ());
    }
 
    /** Calculates the dot product of two vectors.
     *
     * @param a the first Vector3D
     * @param b the second Vector3D
     * @return a Vector3D instance representing the vector product of the two vectors.
     */
 
    public static Vector3D crossVector(Vector3D a, Vector3D b){
    return new Vector3D(a.getY()*b.getZ()-a.getZ()*b.getY(),
                a.getZ()*b.getX()-a.getX()*b.getZ(),
                a.getX()*b.getY()-a.getY()*b.getX());
    }
     
     /** Checks if vectors are equal.
     *
     * @param a the first Vector3D
     * @param b the second Vector3D
     * @return a boolean showing if equality is satisfied.
     */
 
    public static boolean equalVector(Vector3D a, Vector3D b){
    //Use modulus of difference to avoid false negatives due to rounding.
    //eps is epsilon (a small number) which will be used as a threshold.
    double eps = a.mag() * 1E-10;
    return (((Math.abs(a.getX()-b.getX())) < eps) &&
        ((Math.abs(a.getY()-b.getY())) < eps) &&
        ((Math.abs(a.getZ()-b.getZ()))) < eps);
    }
}
