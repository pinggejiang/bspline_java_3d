package Bspline_Java_3D;
import java.io.IOException;

/*
 *
 * @author pingge
 * May 20th, 2014
 * File usage: main B-spline calculation
 *
 */

public class BSPLINE {
    
    public void bspline_xform_initialize(BSPLINE_Xform bxf, float[] img_origin,
            float[] img_spacing, int[] img_dim, int[] roi_offset, int[] roi_dim, 
            int[] vox_per_rgn)
    {
        int d;
        int i, j, k, p;
        int tx, ty, tz;
        float[] A;
        float[] B;
        float[] C;
        System.out.printf ("\nInitializing...\n");
        for (d = 0; d < 3; d++) 
        {
            /* copy input parameters over */
            bxf.img_origin[d] = img_origin[d];
            bxf.img_spacing[d] = img_spacing[d];
            bxf.img_dim[d] = img_dim[d];
            bxf.roi_offset[d] = roi_offset[d];
            bxf.roi_dim[d] = roi_dim[d];
            bxf.vox_per_rgn[d] = vox_per_rgn[d];
            
            /* grid spacing is in mm */
            bxf.grid_spac[d] = bxf.vox_per_rgn[d] * bxf.img_spacing[d];

            /* rdims is the number of regions */
            bxf.rdims[d] = 1 + (bxf.roi_dim[d] - 1) / bxf.vox_per_rgn[d];

            /* cdims is the number of control points */
            bxf.cdims[d] = 3 + bxf.rdims[d];
        }

            /* total number of control points & coefficients */
            bxf.num_knots = bxf.cdims[0] * bxf.cdims[1] * bxf.cdims[2];
            bxf.num_coeff = bxf.cdims[0] * bxf.cdims[1] * bxf.cdims[2] * 3;
            A = new float[bxf.vox_per_rgn[0] * 4];
            B = new float[bxf.vox_per_rgn[1] * 4];
            C = new float[bxf.vox_per_rgn[2] * 4];
            bxf.coeff = new double[bxf.num_coeff];
            bxf.q_lut = new float[bxf.vox_per_rgn[0] * bxf.vox_per_rgn[1] * bxf.vox_per_rgn[2] * 64];

        for (i = 0; i < bxf.vox_per_rgn[0]; i++) 
        {
            float ii = ((float) i) / bxf.vox_per_rgn[0];
            float t3 = ii*ii*ii;
            float t2 = ii*ii;
            float t1 = ii;
            A[i*4+0] = (float) ((1.0/6.0) * (- 1.0 * t3 + 3.0 * t2 - 3.0 * t1 + 1.0));
            A[i*4+1] = (float) ((1.0/6.0) * (+ 3.0 * t3 - 6.0 * t2            + 4.0));
            A[i*4+2] = (float) ((1.0/6.0) * (- 3.0 * t3 + 3.0 * t2 + 3.0 * t1 + 1.0));
            A[i*4+3] = (float) ((1.0/6.0) * (+ 1.0 * t3));
        }

        for (j = 0; j < bxf.vox_per_rgn[1]; j++) 
        {
            float jj = ((float) j) / bxf.vox_per_rgn[1];
            float t3 = jj*jj*jj;
            float t2 = jj*jj;
            float t1 = jj;
            B[j*4+0] = (float)((1.0/6.0) * (- 1.0 * t3 + 3.0 * t2 - 3.0 * t1 + 1.0));
            B[j*4+1] = (float)((1.0/6.0) * (+ 3.0 * t3 - 6.0 * t2            + 4.0));
            B[j*4+2] = (float)((1.0/6.0) * (- 3.0 * t3 + 3.0 * t2 + 3.0 * t1 + 1.0));
            B[j*4+3] = (float)((1.0/6.0) * (+ 1.0 * t3));
        }

        for (k = 0; k < bxf.vox_per_rgn[2]; k++) 
        {
            float kk = ((float) k) / bxf.vox_per_rgn[2];
            float t3 = kk*kk*kk;
            float t2 = kk*kk;
            float t1 = kk;
            C[k*4+0] = (float)((1.0/6.0) * (- 1.0 * t3 + 3.0 * t2 - 3.0 * t1 + 1.0));
            C[k*4+1] = (float)((1.0/6.0) * (+ 3.0 * t3 - 6.0 * t2            + 4.0));
            C[k*4+2] = (float)((1.0/6.0) * (- 3.0 * t3 + 3.0 * t2 + 3.0 * t1 + 1.0));
            C[k*4+3] = (float)((1.0/6.0) * (+ 1.0 * t3));
        }

        p = 0;
        for (k = 0; k < bxf.vox_per_rgn[2]; k++) 
        {
            for (j = 0; j < bxf.vox_per_rgn[1]; j++) 
            {
                for (i = 0; i < bxf.vox_per_rgn[0]; i++) 
                {
                    for (tz = 0; tz < 4; tz++) {
                        for (ty = 0; ty < 4; ty++) {
                            for (tx = 0; tx < 4; tx++) {
                                bxf.q_lut[p++] = A[i*4+tx] * B[j*4+ty] * C[k*4+tz];
                            }
                        }
                    }
                }
            }
        }
        p = 0;
        bxf.c_lut = new int[bxf.rdims[0] * bxf.rdims[1] * bxf.rdims[2] * 64];
        for (k = 0; k < bxf.rdims[2]; k++) {
            for (j = 0; j < bxf.rdims[1]; j++) {
                for (i = 0; i < bxf.rdims[0]; i++) {
                    for (tz = 0; tz < 4; tz++) {
                        for (ty = 0; ty < 4; ty++) {
                            for (tx = 0; tx < 4; tx++) {
                                bxf.c_lut[p++] = 
                                    + (k + tz) * bxf.cdims[0] * bxf.cdims[1]
                                    + (j + ty) * bxf.cdims[0] 
                                    + (i + tx);
                            }
                        }
                    }
                }
            }
        }
        System.out.printf ("            # of Tiles:  %d\t[%d x %d x %d]\n", 
        bxf.rdims[0] * bxf.rdims[1] * bxf.rdims[2],
        bxf.rdims[0], bxf.rdims[1], bxf.rdims[2]);

        System.out.printf ("       Voxels per Tile:  %d\t[%d x %d x %d]\n", 
        bxf.vox_per_rgn[0] * bxf.vox_per_rgn[1] * bxf.vox_per_rgn[2],
        bxf.vox_per_rgn[0], bxf.vox_per_rgn[1], bxf.vox_per_rgn[2]);

        System.out.printf ("   # of Control Points:  %d\t[%d x %d x %d]\n\n", 
        bxf.cdims[0] * bxf.cdims[1] * bxf.cdims[2],
        bxf.cdims[0], bxf.cdims[1], bxf.cdims[2]);

    }
    public void bspline_run_optimization(BSPLINE_Xform bxf, BSPLINE_Parms parms,
            Volume fixed, Volume moving, Volume moving_grad) throws IOException
    {
        BSPLINE_State bst = new BSPLINE_State();
        bspline_optimize(bxf, bst, parms, fixed, moving, moving_grad);
    }
    
    
    public void bspline_optimize(BSPLINE_Xform bxf, BSPLINE_State bst,
            BSPLINE_Parms parms, Volume fixed, Volume moving, Volume moving_grad) throws IOException
    {
        Bspline_optimize_data bod = new Bspline_optimize_data();
        bod.bxf = bxf;
        bod.bst = bst;
        bod.parms = parms;
        bod.fixed = fixed;
        bod.moving = moving;
        bod.moving_grad = moving_grad;
        int n = bod.bxf.num_coeff;
        int m = 6;
        int[] iprint = new int[2];
        iprint[ 1 -1] = 1;
        iprint[ 2 -1] = 0;
        boolean diagco = false;
        double diag[] = new double[n];
        double eps = 1.0e-5;
        double xtol = 1.0e-16;
        int icall = 0;
        int[] iflag = new int[1];
        iflag[0] = 0;
        double old_ssd = 0;
        BSPLINE_SCORE_CAL Bspline_score_cal = new BSPLINE_SCORE_CAL();
        do
        { 
            old_ssd = bod.bst.ssd.score;
            Bspline_score_cal.bspline_score_cal(bod.parms, bod.bst, bod.bxf, bod.fixed, bod.moving, bod.moving_grad);  
            try
            {
                LBFGS.lbfgs(n, m, bod.bxf.coeff, bod.bst.ssd.score, bod.bst.ssd.grad, diagco, diag, iprint, eps, xtol, iflag);
                System.out.printf("Complete one Iteration\n" + "SSD: %f\n", bod.bst.ssd.score);
            }
            catch (LBFGS.ExceptionWithIflag e)
            {
                    System.err.println( "Sdrive: lbfgs failed.\n"+e );
                    return;
            }
            
            icall += 1;
        }
        while (iflag[0] != 0 && icall <= 200 && Math.abs(bod.bst.ssd.score - old_ssd) >= 5);
    }

    public void bspline_interpolate_vf (Volume interp, BSPLINE_Xform bxf)
    {
        int i, j, k, v;
        int[] p = new int[3];
        int[] q = new int[3];
        float[] out = new float[3];
        int qidx;

        for (k = 0; k < bxf.roi_dim[2]; k++) 
        {
            p[2] = k / bxf.vox_per_rgn[2];
            q[2] = k % bxf.vox_per_rgn[2];
            for (j = 0; j < bxf.roi_dim[1]; j++) 
            {
                p[1] = j / bxf.vox_per_rgn[1];
                q[1] = j % bxf.vox_per_rgn[1];
                for (i = 0; i < bxf.roi_dim[0]; i++) 
                {
                    p[0] = i / bxf.vox_per_rgn[0];
                    q[0] = i % bxf.vox_per_rgn[0];
                    qidx = (((q[2] * bxf.vox_per_rgn[1] + q[1]) * bxf.vox_per_rgn[0]) + q[0]);
                    v = (k+bxf.roi_offset[2]) * interp.dim[0] * interp.dim[1] + (j+bxf.roi_offset[1]) * interp.dim[0] 
                        + (i+bxf.roi_offset[0]);
                    out = bspline_interp_pix (out, bxf, p, qidx);
                    interp.img[3*v] = out[0];
                    interp.img[3*v + 1] = out[1];
                    interp.img[3*v + 2] = out[2];
                }
            }
        }
    }
    

    public float[] bspline_interp_pix (float[] out, BSPLINE_Xform bxf, int[] p, int qidx)
    {
        int i, j, k, m;
        int cidx;
        float[] q_lut = new float[64];
        for (i = 0; i < 64; i++)
        {
            q_lut[i] = bxf.q_lut[qidx*64+i];
        }
        out[0] = out[1] = out[2] = 0;
        m = 0;
        for (k = 0; k < 4; k++) 
        {
            for (j = 0; j < 4; j++) 
            {
                for (i = 0; i < 4; i++) 
                {
                    cidx = (p[2] + k) * bxf.cdims[1] * bxf.cdims[0]
                        + (p[1] + j) * bxf.cdims[0]
                        + (p[0] + i);
                    cidx = cidx * 3;
                    out[0] += q_lut[m] * bxf.coeff[cidx+0];
                    out[1] += q_lut[m] * bxf.coeff[cidx+1];
                    out[2] += q_lut[m] * bxf.coeff[cidx+2];
                    m ++;
                }
            }
        }
        return out;
    }
    
}
