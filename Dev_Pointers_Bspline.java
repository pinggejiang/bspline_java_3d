package Bspline_Java_3D;

/*
 *
 * @author pingge
 * May 20th, 2014
 * File usage: class definition
 *
 */

public class Dev_Pointers_Bspline {

    float[] fixed_image;     
    float[] moving_image;   
    float[] moving_grad;     

    float[] coeff;        
    float[] score;       

    float[] dc_dv;       
    float[] dc_dv_x;       
    float[] dc_dv_y;   
    float[] dc_dv_z;        

    float[] cond_x;       
    float[] cond_y;          
    float[] cond_z;         

    float[] grad;        
    float[] dc_dp_x;
    float[] dc_dp_y;
    float[] dc_dp_z;
    float[] grad_temp;

    int[] LUT_Knot;
    int[] LUT_NumTiles;
    int[] LUT_Offsets;
    float[] LUT_Bspline_x;
    float[] LUT_Bspline_y;
    float[] LUT_Bspline_z;
    float[] skipped;    

}
