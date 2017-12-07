package Bspline_Java_3D;

/*
 *
 * @author pingge
 * May 20th, 2014
 * File usage: class definition
 *
 */

public class BSPLINE_Xform {
    float[] img_origin = new float[3];         /* Image origin (in mm) */
    float[] img_spacing = new float[3];        /* Image spacing (in mm) */
    int[] img_dim = new int[3];              /* Image size (in vox) */
    int[] roi_offset = new int[3];           /* Position of first vox in ROI (in vox) */
    int[] roi_dim = new int[3];              /* Dimension of ROI (in vox) */
    int[] vox_per_rgn = new int[3];          /* Knot spacing (in vox) */
    float[] grid_spac = new float[3];          /* Knot spacing (in mm) */
    int[] rdims = new int[3];                /* # of regions in (x,y,z) */
    int[] cdims = new int[3];                /* # of knots in (x,y,z) */
    int num_knots;               /* Total number of knots (= product(cdims)) */
    int num_coeff;               /* Total number of coefficents (= product(cdims) * 3) */
    double[] coeff;                /* Coefficients.  Vector directions interleaved. */
    int[] cidx_lut;               /* Lookup volume for region number */
    int[] c_lut;                  /* Lookup table for control point indices */
    int[] qidx_lut;               /* Lookup volume for region offset */
    float[] q_lut;                /* Lookup table for influence multipliers */
}
