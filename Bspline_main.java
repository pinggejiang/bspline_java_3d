package Bspline_Java_3D;
import java.io.IOException;

/*
 *
 * @author pingge
 * May 20th, 2014
 * File usage: main class
 *
 */

//"Usage: bspline [options] fixed_image moving_image\n"
//"Options:\n"
/* " -m iterations              Maximum iterations (default is 10)\n"
    " -s \"i j k\"              Integer knot spacing (voxels)\n"
    " -V outfile                Output vector field\n"
    " -O outfile                Output warped image\n"
    " -fixed 
    " -moving
*/
public class Bspline_main {
    public static void main(String[] args) throws IOException {
        /*global parameter init*/
        BSPLINE_Options options = new BSPLINE_Options();
        BSPLINE_Parms parms = new BSPLINE_Parms();
        BSPLINE_Xform bxf = new BSPLINE_Xform();
        Volume moving;
        Volume fixed;
        Volume moving_grad = new Volume();
        Volume vector_field = new Volume();
        Volume moving_warped;
        int[] roi_offset = new int[3]; // used only if roi is not the whole image
        BSPLINE_Opts opts = new BSPLINE_Opts();
        /* reading input parameters*/
        opts.bspline_opts_parse_args(options, args);
        mha_oper mha = new mha_oper();
        fixed = mha.Read_mha(options.fixed_fn);
        System.out.println("Done with reading fixed image");
        moving = mha.Read_mha(options.moving_fn);
        System.out.println("Done with reading moving image");
        /*compute moving image gradient*/
        System.out.printf ("Computing Moving Image Gradient...\n");
        Volume_proc vol_proc = new Volume_proc();
        moving_grad = vol_proc.volume_make_gradient (moving);
        System.out.printf ("Allocating Lookup Tables...\n");
        /*init B-spline transform parameters*/
        BSPLINE bs = new BSPLINE();
        bs.bspline_xform_initialize(bxf, fixed.offset, fixed.pix_spacing,
                fixed.dim, roi_offset, fixed.dim, options.vox_per_rgn);
        /* Run the optimization */
        System.out.printf ("Running optimization.\n");
        bs.bspline_run_optimization (bxf, parms, fixed, moving, moving_grad);
        System.out.printf ("Done running optimization.\n");
        /* Create vector field from bspline coefficients and save */
        if (options.output_vf_fn != null || options.output_warped_fn!= null)
        {
            System.out.printf ("Creating vector field.\n");
            vector_field = vol_proc.volume_create (fixed.dim, fixed.offset, fixed.pix_spacing, "PT_VF_FLOAT_INTERLEAVED",
                    fixed.direction_cosines, 0);
            bs.bspline_interpolate_vf (vector_field, bxf);
        }
        /* Create warped output image and save */
        if (options.output_warped_fn != null) 
        {
            System.out.printf ("Warping image.\n");
            moving_warped = vol_proc.vf_warp (null, moving, vector_field);
            if (moving_warped != null) 
            {
                System.out.printf ("Writing warped image.\n");
                mha.write_mha (options.output_warped_fn, moving_warped);
            } else 
            {
                System.err.printf ("Sorry, couldn't create warped image.\n");
            }
        }
    }
}
