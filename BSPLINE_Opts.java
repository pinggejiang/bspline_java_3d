package Bspline_Java_3D;

/*
 *
 * @author pingge
 * May 20th, 2014
 * File usage: receive input arguments
 *
 */
public class BSPLINE_Opts {
    //"Usage: bspline [options] \n"
    
/*  " -m iterations              Maximum iterations (default is 10)\n"
    " -s \"i j k\"                 Integer knot spacing (voxels)\n"
    " -V outfile                 Output vector field\n"
    " -O outfile                 Output warped image\n"
    */
    public void bspline_opts_parse_args (BSPLINE_Options options, String[] args)
    {
        int d, i, j;
        int[] rc;
        BSPLINE_Parms parms = new BSPLINE_Parms();
        parms = options.parms;
        for (d = 0; d <3; d++)
        {
            options.vox_per_rgn[d] = 15;
        }
        bspline_parms_set_default(parms);
        try
        {
            for (i = 0; i < args.length; i++)
            {
                switch (args[i]) {
                    case "-s":
                        options.vox_per_rgn[0] = Integer.parseInt(args[i+1]);
                        options.vox_per_rgn[1] = Integer.parseInt(args[i+2]);
                        options.vox_per_rgn[2] = Integer.parseInt(args[i+3]);
                        if (args[i+2].substring(0).equals("-"))
                        {
                            options.vox_per_rgn[0] = Integer.parseInt(args[i+1]);
                            options.vox_per_rgn[1] = Integer.parseInt(args[i+1]);
                            options.vox_per_rgn[2] = Integer.parseInt(args[i+1]);
                        }
                        else if (args[i+3].substring(0).equals("-"))
                        {
                            System.out.println("option %s requires three arguments");
                        }   break;
                    case "-O":
                        options.output_warped_fn = args[i+1];
                        break;
                    case "-v":
                        options.output_vf_fn = args[i+1];
                        break;
                    case "-fixed":
                        options.fixed_fn = args[i+1];
                        break;
                    case "-moving":
                        options.moving_fn = args[i+1];
                        break;
                }
            }
            if (args.length <2)
            {
                String str = "\" -m iterations             Maximum iterations (default is 10)\\n\"\n" +
                        "    \" -s \\\"i j k\\\"           Integer knot spacing (voxels)\\n\"\n" +
                        "    \" -v outfile                 Output vector field\\n\"\n" +
                        "    \" -O outfile                 Output warped image\\n\"\n" +
                        "    \" -fixed image               Input fixed image\\n\"\n" +
                        "    \" -moving image              Input moving image\\n\"";
                System.out.println (str);
                System.exit(1);
            }
            System.out.println("Preparing to register the following images: \n");
            System.out.printf("    *   Fixed: %s\n", options.fixed_fn);
            System.out.printf("    *   Moving: %s\n", options.moving_fn);
            System.out.println("\n");
        }
        catch(Exception e)
        {
            String str = "\" -m iterations              Maximum iterations (default is 10)\\n\"\n" +
                        "    \" -s \\\"i j k\\\"        Integer knot spacing (voxels)\\n\"\n" +
                        "    \" -v outfile              Output vector field\\n\"\n" +
                        "    \" -O outfile              Output warped image\\n\"\n" +
                        "    \" -fixed image            Input fixed image\\n\"\n" +
                        "    \" -moving image           Input moving image\\n\"";
                System.err.println (str);
                System.exit(1);
        }
        
    }
    
    public void bspline_parms_set_default(BSPLINE_Parms parms)
    {
        parms.max_its = 10;
        parms.convergence_tol = 0.1;
        parms.convergence_tol_its = 4;
    }

}




