package Bspline_Java_3D;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.lang.Math;

/*
 *
 * @author pingge
 * May 20th, 2014
 * File usage: B-spline score calculation and belonging sub_functions
 *
 */

public class BSPLINE_SCORE_CAL {
    public void bspline_score_cal(BSPLINE_Parms parms, BSPLINE_State bst, BSPLINE_Xform bxf,
            Volume fixed, Volume moving, Volume moving_grad) throws IOException 
    {
        BSPLINE_Score ssd = bst.ssd;
        double score_tile;
        int num_vox;

        float[] f_img = fixed.img;
        float[] m_img = moving.img;
        float[] m_grad = moving_grad.img;

        int idx_tile;
        int num_tiles = bxf.rdims[0] * bxf.rdims[1] * bxf.rdims[2];

        
        float[] cond_x = new float[64*bxf.num_knots];
        float[] cond_y = new float[64*bxf.num_knots];
        float[] cond_z = new float[64*bxf.num_knots];
        
        int idx_knot;
        int idx_set;
        int i;
        float[] q_lut;
        ssd.score = 0;
        num_vox = 0;
        score_tile = 0;
        ssd.grad = new double[bxf.num_coeff];
        for (idx_tile = 0; idx_tile < num_tiles; idx_tile++) {
            int rc;
            int set_num;
            int[] crds_tile =new int[3];
            int[] crds_local = new int[3];
            int idx_local;

            float[] phys_fixed = new float[3];
            int[] crds_fixed = new int[3];
            int idx_fixed;

            float[] dxyz = new float[3];

            float[] phys_moving = new float[3];
            float[] crds_moving = new float[3];
            int[] crds_moving_floor = new int[3];
            int[] crds_moving_round = new int[3];
            int idx_moving_floor;
            int idx_moving_round;

            float[] li_1 = new float[3];
            float[] li_2 = new float[3];
            float m_val = 0, diff;

            float[] dc_dv = new float[3];

            float[] sets_x = new float[64];
            float[] sets_y = new float[64];
            float[] sets_z = new float[64];

            int[] k_lut = new int[64];

            // Get tile coordinates from index
             COORDS_FROM_INDEX (crds_tile, idx_tile, bxf.rdims);

            // Serial through voxels in tile
            for (crds_local[2] = 0; crds_local[2] < bxf.vox_per_rgn[2]; crds_local[2]++) {
                for (crds_local[1] = 0; crds_local[1] < bxf.vox_per_rgn[1]; crds_local[1]++) {
                    for (crds_local[0] = 0; crds_local[0] < bxf.vox_per_rgn[0]; crds_local[0]++) {
                        q_lut = new float[64];
                        crds_fixed[0] = bxf.roi_offset[0] + bxf.vox_per_rgn[0] * crds_tile[0] + crds_local[0];
                        crds_fixed[1] = bxf.roi_offset[1] + bxf.vox_per_rgn[1] * crds_tile[1] + crds_local[1];
                        crds_fixed[2] = bxf.roi_offset[2] + bxf.vox_per_rgn[2] * crds_tile[2] + crds_local[2];

                        // Make sure we are inside the image volume
                        if (crds_fixed[0] >= bxf.roi_offset[0] + bxf.roi_dim[0]) 
                    	{
                        	continue;
                    	}
                        if (crds_fixed[1] >= bxf.roi_offset[1] + bxf.roi_dim[1]) 
                    	{
                        	continue;
                        }
                        if (crds_fixed[2] >= bxf.roi_offset[2] + bxf.roi_dim[2]) 
                    	{
                        	continue;
                    	}

                        // Compute physical coordinates of fixed image voxel
                        phys_fixed[0] = bxf.img_origin[0] + bxf.img_spacing[0] * crds_fixed[0];
                        phys_fixed[1] = bxf.img_origin[1] + bxf.img_spacing[1] * crds_fixed[1];
                        phys_fixed[2] = bxf.img_origin[2] + bxf.img_spacing[2] * crds_fixed[2];

                        // Construct the local index within the tile
                        idx_local = (((crds_local[2] * bxf.vox_per_rgn[1] + crds_local[1]) * bxf.vox_per_rgn[0]) + crds_local[0]);
                        // Construct the image volume index
                        idx_fixed = (((crds_fixed[2] * fixed.dim[1] + crds_fixed[1]) * fixed.dim[0]) + crds_fixed[0]);
                        // Calc. deformation vector (dxyz) for voxel
                        bspline_interp_pix_b (dxyz, bxf, idx_tile, idx_local);
                        
                        /* some variable name definitions:
                        *crds_tile: tile index within the volume
                        *crds_local: voxel coordinates within the tile (same for each tile)
                        *crds_fixed: voxel coordinates within the volume
                        *phys_fixed: voxel coordinates in the world
                        *idx_local: voxel index within the tile
                        *idx_fixed: voxel index within the volume
                         */
                        
                        // Calc. moving image coordinate from the deformation vector
                        rc = bspline_find_correspondence (phys_moving, crds_moving, phys_fixed, dxyz, moving);

                        // Return code is 0 if voxel is pushed outside of moving image
                        if (rc == 0 ) 
                    	{
                        	continue;
                    	}

                        // Compute linear interpolation fractions
                        CLAMP_LINEAR_INTERPOLATE_3D (crds_moving, crds_moving_floor, crds_moving_round, li_1, li_2, moving);

                        // Find linear indices for moving image
                        idx_moving_floor = (((crds_moving_floor[2] * moving.dim[1] + crds_moving_floor[1]) * moving.dim[0]) + crds_moving_floor[0]);
                        idx_moving_round = (((crds_moving_round[2] * moving.dim[1] + crds_moving_round[1]) * moving.dim[0]) + crds_moving_round[0]);

                        // Calc. moving voxel intensity via linear interpolation
                        m_val = BSPLINE_LI_VALUE (m_val, li_1[0], li_2[0], li_1[1], li_2[1], li_1[2], li_2[2], idx_moving_floor, m_img, moving);

                        // Compute intensity difference
                        diff = m_val - f_img[idx_fixed];

                        // Store the score!
                        score_tile += diff * diff;
                        num_vox++;

                        // Compute dc_dv
                        dc_dv[0] = diff * m_grad[3 * idx_moving_round + 0];
                        dc_dv[1] = diff * m_grad[3 * idx_moving_round + 1];
                        dc_dv[2] = diff * m_grad[3 * idx_moving_round + 2];

                        // Initialize q_lut
                        for (i = 0; i < 64; i++)
                        {
                        	q_lut[i] = bxf.q_lut[64 * idx_local +i];
                        }
                        // Condense dc_dv @ current voxel index
                        for (set_num = 0; set_num < 64; set_num++) {
                            sets_x[set_num] += dc_dv[0] * q_lut[set_num];
                            sets_y[set_num] += dc_dv[1] * q_lut[set_num];
                            sets_z[set_num] += dc_dv[2] * q_lut[set_num];
                        }
                    }
                }
            }//done with loop through voxels within tiles
            find_knots(k_lut, idx_tile, bxf.rdims, bxf.cdims);

            for (set_num = 0; set_num < 64; set_num++) {
                int knot_num = k_lut[set_num];
                cond_x[ (64*knot_num) + (63 - set_num) ] = sets_x[set_num];
                cond_y[ (64*knot_num) + (63 - set_num) ] = sets_y[set_num];
                cond_z[ (64*knot_num) + (63 - set_num) ] = sets_z[set_num];
            }
        }// done with loop through tiles

        for (idx_knot = 0; idx_knot < (bxf.cdims[0] * bxf.cdims[1] * bxf.cdims[2]); idx_knot++) {
            for(idx_set = 0; idx_set < 64; idx_set++) {
                ssd.grad[3*idx_knot + 0] += cond_x[64*idx_knot + idx_set];
                ssd.grad[3*idx_knot + 1] += cond_y[64*idx_knot + idx_set];
                ssd.grad[3*idx_knot + 2] += cond_z[64*idx_knot + idx_set];
            }
        }
        ssd.score = score_tile / num_vox;

        for (i = 0; i < bxf.num_coeff; i++) {
            ssd.grad[i] = 2 * ssd.grad[i] / num_vox;
        }
    }
    
    public void COORDS_FROM_INDEX (int[] ijk, int idx, int[] dim)
    {
        ijk[2] = idx / (dim[0] * dim[1]);
        ijk[1] = (idx-(ijk[2] * dim[0] * dim[1])) / dim[0];
        ijk[0] = idx - ijk[2] * dim[0] * dim[1] - (ijk[1] * dim[0]);
    }
    
    public void bspline_interp_pix_b(float[] out, BSPLINE_Xform bxf, int pidx, int qidx)
    {
        int i, j, k, m;
        int cidx;
        float[] q_lut = new float[64];
        for (i = 0; i< 64; i++)
        {
        	q_lut[i] = bxf.q_lut[qidx * 64 + i];
        }
    	int[] c_lut =new int[64];
        for (i = 0; i <64; i++)
    	{
        	c_lut[i] = bxf.c_lut[pidx * 64 + i];
    	}

        out[0] = out[1] = out[2] = 0;
        m = 0;
        for (k = 0; k < 4; k++) 
        {
            for (j = 0; j < 4; j++) 
            {
                for (i = 0; i < 4; i++) 
                {
                    cidx = 3 * c_lut[m];
                    out[0] += q_lut[m] * bxf.coeff[cidx+0];
                    out[1] += q_lut[m] * bxf.coeff[cidx+1];
                    out[2] += q_lut[m] * bxf.coeff[cidx+2];
                    m ++;
                }
            }
        }
    }
    
   public int bspline_find_correspondence 
   (
    float[] mxyz,      /* xyz coordinates in moving image (mm) */
    float[] mijk,      /* ijk indices in moving image (vox) */
    float[] fxyz,      /* xyz coordinates in fixed image (mm) */
    float[] dxyz,      /* displacement from fixed to moving (mm) */
    Volume moving      /* moving image */
    )
   {
       mxyz[0] = fxyz[0] + dxyz[0];
       mijk[0] = (mxyz[0] - moving.offset[0]) / moving.pix_spacing[0];
       mxyz[1] = fxyz[1] + dxyz[1];
       mijk[1] = (mxyz[1] - moving.offset[1]) / moving.pix_spacing[1];
       mxyz[2] = fxyz[2] + dxyz[2];
       mijk[2] = (mxyz[2] - moving.offset[2]) / moving.pix_spacing[2];
       if (mijk[0] < -0.5 || mijk[0] > moving.dim[0] - 0.5 || mijk[1] < -0.5 || mijk[1] > 
       moving.dim[1] - 0.5 || mijk[2] < -0.5 || mijk[2] > moving.dim[2] - 0.5) 
	   {
    	   return 0;
	   }
       else
       {
    	   return 1;
       }
   }
   
   public void CLAMP_LINEAR_INTERPOLATE_3D(float[] mijk, int[] mijk_f, int[] mijk_r, float[] li_frac_1, 
           float[] li_frac_2, Volume moving)
   {
       for (int i = 0; i < 3; i++)
       {
    	   float maff = (float)Math.floor(mijk[i]);
    	   mijk_f[i] = (int)maff;
    	   if (mijk[i] >= 0)
    	   {
    		   mijk_r[i] = (int)((mijk[i]) + 0.5);
    	   }else
    	   {
    		   mijk_r[i] = (int)(-(-(mijk[i]) + 0.5));
    	   }
    	   li_frac_2[i] = mijk[i] - maff;
    	   if (mijk_f[i] < 0)
    	   {
    		   mijk_f[i] = 0;
    		   mijk_r[i] = 0;
    		   li_frac_2[i] = 0.0f;
    	   }else if (mijk_f[i] >= moving.dim[i] - 1)
    	   {
    		   mijk_f[i] = moving.dim[i] - 2;
    		   mijk_r[i] = moving.dim[i] - 1;
    		   li_frac_2[i] = 1.0f;
    	   }
    	   li_frac_1[i] = 1.0f - li_frac_2[i];
       }
   }
   
   public void clamp_linear_interpolate (
            float ma,          /*  Unrounded pixel coordinate */
            int dmax,          /* Maximum coordinate */
            int maf,           /* floor pixel*/
            int mar,           /* round pixel */
            float fa1,         /* Fraction for lower index voxel */
            float fa2          /* Fraction for upper index voxel */
    )
    {
        float maff = (float)Math.floor(ma);
        maf = (int) maff;
        if (ma >= 0)
        {
            mar = (int)((ma)+0.5);
        }
        else
        {
            mar = (int)(-(-(ma)+0.5));
        }
        fa2 = ma - maff;

        if (maf < 0) {
            maf = 0;
            mar = 0;
            fa2 = 0.0f;
        } else if (maf >= dmax) {
            maf = dmax - 1;
            mar = dmax;
            fa2 = 1.0f;
        }

        fa1 = 1.0f - fa2;
    }
    
   public float BSPLINE_LI_VALUE(float m_val, float fx1, float fx2, float fy1, float fy2, float fz1, float fz2, int mvf,
           float[] m_img, Volume moving)
   {
    						   
	float m_x1y1z1, m_x2y1z1, m_x1y2z1, m_x2y2z1;		   
	float m_x1y1z2, m_x2y1z2, m_x1y2z2, m_x2y2z2;		   
								   
	m_x1y1z1 = fx1 * fy1 * fz1 * m_img[mvf];		   
	m_x2y1z1 = fx2 * fy1 * fz1 * m_img[mvf+1];		   
	m_x1y2z1 = fx1 * fy2 * fz1 * m_img[mvf+moving.dim[0]];		
	m_x2y2z1 = fx2 * fy2 * fz1 * m_img[mvf+moving.dim[0]+1];	
	m_x1y1z2 = fx1 * fy1 * fz2 * m_img[mvf+moving.dim[1]*moving.dim[0]]; 
	m_x2y1z2 = fx2 * fy1 * fz2 * m_img[mvf+moving.dim[1]*moving.dim[0]+1]; 
	m_x1y2z2 = fx1 * fy2 * fz2 * m_img[mvf+moving.dim[1]*moving.dim[0]+moving.dim[0]]; 
	m_x2y2z2 = fx2 * fy2 * fz2 * m_img[mvf+moving.dim[1]*moving.dim[0]+moving.dim[0]+1]; 
	m_val = m_x1y1z1 + m_x2y1z1 + m_x1y2z1 + m_x2y2z1 + m_x1y1z2 + m_x2y1z2 + m_x1y2z2 + m_x2y2z2;	
	return m_val;
   }
   
   public void find_knots(int[] knots, int tile_num, int[] rdims, int[] cdims)
   {
        int[] tile_loc = new int[3];
        int i, j, k;
        int idx = 0;
        int num_tiles_x = cdims[0] - 3;
        int num_tiles_y = cdims[1] - 3;
        int num_tiles_z = cdims[2] - 3;

        tile_loc[0] = tile_num % num_tiles_x;
        tile_loc[1] = ((tile_num - tile_loc[0]) / num_tiles_x) % num_tiles_y;
        tile_loc[2] = ((((tile_num - tile_loc[0]) / num_tiles_x) / num_tiles_y) % num_tiles_z);
        
        tile_loc[0]++;
        tile_loc[1]++;
        tile_loc[2]++;

        for (k = -1; k < 3; k++) {
            for (j = -1; j < 3; j++) {
                for (i = -1; i < 3; i++) {
                    knots[idx++] = (cdims[0]*cdims[1]*(tile_loc[2]+k)) + (cdims[0]*(tile_loc[1]+j)) + (tile_loc[0]+i);
                }
            }
        }
   }

}
