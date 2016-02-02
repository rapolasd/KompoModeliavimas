1
2
3
4
5
6
7
8
9
10
11
12
13
14
15
16
17
18
19
20
21
22
23
24
25
26
27
28
29
30
31
32
33
34
35
36
37
38
39
40
41
42
43
44
45
46
47
48
49
50
51
52
53
54
55
56
57
58
59
60
61
62
63
64
65
66
67
68
69
70
71
72
73
74
75
76
77
78
79
80
81
82
83
84
85
86
87
88
89
90
91
92
93
94
95
96
97
98
99
100
101
102
103
104
105
106
107
108
109
110
111
112
113
114
115
116
117
118
119
120
121
122
123
124
125
126
127
128
129
130
131
132
133
134
135
136
137
138
139
140
141
142
143
144
145
146
147
148
149
150
151
152
153
154
155
156
157
158
159
160
161
162
163
164
165
166
167
168
169
170
171
172
173
174
175
176
177
178
179
180
181
182
183
184
185
186
187
188
189
190
191
192
193
194
195
196
197
198
199
200
201
202
203
204
205
206
207
208
209
210
211
212
213
214
215
216
217
218
219
220
221
222
223
224
225
226
227
228
229
230
231
232
233
234
235
236
237
238
239
240
241
242
243
244
245
246
247
248
249
250
251
252
253
254
255
256
257
258
259
260
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
     
     
    /*
     * Static methods
     *
     *
     */
 
    /** Adds two vectors together.
     *
     * @param a the first Vector3D
     * @param b the second Vector3D
     * @return a vector3D representing the sum of the two vectors.
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
     * @return a Vector3D representing first Vector3D subtracted by a second Vector3D.
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
     * @return a Vector3D representing the vector product of the two vectors.
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
