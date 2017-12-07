package Bspline_Java_3D;

/*
 *
 * @author pingge
 * May 20th, 2014
 * File usage: class definition
 *
 */
public class Bspline_optimize_data {
    BSPLINE_Xform bxf = new BSPLINE_Xform();
    BSPLINE_State bst = new BSPLINE_State();
    BSPLINE_Parms parms = new BSPLINE_Parms();
    Volume fixed = new Volume();
    Volume moving = new Volume();
    Volume moving_grad = new Volume();
}
