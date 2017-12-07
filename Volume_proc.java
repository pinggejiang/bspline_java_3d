package Bspline_Java_3D;
import java.util.Arrays;

/*
 *
 * @author pingge
 * May 20th, 2014
 * File usage: volume class macros
 * functions(methods) included:
 * volume_make_gradient(main)
 * volume_create
 * volume_calc_grad(sub)
 * volume_index(get index of a pixel)
 * vf_warp(warp volume based on estimated transformation)
 * round_int
 *
 */

public class Volume_proc {
    public Volume volume_make_gradient (Volume ref)
    {
        Volume grad = volume_create(ref.dim, ref.offset, ref.pix_spacing,
                                    "PT_VF_FLOAT_INTERLEAVED", ref.direction_cosines, 0);
        volume_calc_grad(grad, ref);
        return grad;
    }
    
    public Volume volume_create(int[] dim, float[] offset, float[] pix_spacing,
        String pix_type, float[] direction_cosines, int min_size)
    {
        Volume vol = new Volume();
        int i;
        for (i = 0; i < 3; i++)
        {
            vol.dim[i] = dim[i];
            vol.offset[i] = offset[i];
            vol.pix_spacing[i] = pix_spacing[i];
        }
        if (direction_cosines !=null)
        {
            vol.direction_cosines = direction_cosines;
        }
        else
        {
            vol.direction_cosines[0] = 1.0f;
            vol.direction_cosines[4] = 1.0f;
            vol.direction_cosines[8] = 1.0f;
        }
        vol.npix = vol.dim[0] * vol.dim[1] * vol.dim[2];
        vol.pix_type = pix_type;
        switch (pix_type) 
        {
            case "PT_UCHAR":
                vol.pix_size = 1;
                vol.img = new float[vol.npix]; 
                break;
            case "PT_SHORT":
                vol.pix_size = 2;
                vol.img = new float[vol.npix]; 
                break;
            case "PT_UINT16":
                vol.pix_size = 2;
                vol.img = new float[vol.npix]; 
                break;
            case "PT_UINT32":
                vol.pix_size = 4;
                vol.img = new float[vol.npix]; 
                break;
            case "PT_FLOAT":
                vol.pix_size = 4;
                vol.img = new float[vol.npix]; 
                break;
            case "PT_VF_FLOAT_INTERLEAVED":
                vol.pix_size = 3 * 4;
                vol.img = new float[vol.npix * 3]; 
                break;
            default:
                System.out.printf ("Unhandled type in volume_create().\n");
                System.exit (-1);
        }
        Arrays.fill(vol.img, 0.0f);
        return vol;
    }
    
    /*volume gradient*/
    public void volume_calc_grad(Volume vout, Volume vref)
    {
        int v;
        int i_p, i, i_n, j_p, j, j_n, k_p, k, k_n; /* p is prev, n is next */
        int gi, gj, gk;
        int idx_p, idx_n;
        float[] out_img;
        float[] ref_img;
        out_img = vout.img;
        ref_img = vref.img;
        v = 0;
        for (k = 0; k < vref.dim[2]; k++)
        {
            k_p = k - 1;
            k_n = k + 1;
            if (k == 0) 
            {
                k_p = 0;
            }
            if (k == vref.dim[2]-1) 
            {
                k_n = vref.dim[2]-1;
            }
            for (j = 0; j < vref.dim[1]; j++) 
            {
                j_p = j - 1;
                j_n = j + 1;
                if (j == 0) 
                {
                    j_p = 0;
                }
                if (j == vref.dim[1]-1) 
                {
                    j_n = vref.dim[1]-1;
                }
                for (i = 0; i < vref.dim[0]; i++, v++) 
                {
                    i_p = i - 1;
                    i_n = i + 1;
                    if (i == 0) 
                    {
                        i_p = 0;
                    }
                    if (i == vref.dim[0]-1) 
                    {
                        i_n = vref.dim[0]-1;
                    }

                    gi = 3 * v + 0;
                    gj = 3 * v + 1;
                    gk = 3 * v + 2;

                    idx_p = volume_index (vref.dim, i_p, j, k);
                    idx_n = volume_index (vref.dim, i_n, j, k);
                    out_img[gi] = (float)( (ref_img[idx_n] - ref_img[idx_p]) / 2.0 / vref.pix_spacing[0]);
                    idx_p = volume_index (vref.dim, i, j_p, k);
                    idx_n = volume_index (vref.dim, i, j_n, k);
                    out_img[gj] = (float)( (ref_img[idx_n] - ref_img[idx_p]) / 2.0 / vref.pix_spacing[1]);
                    idx_p = volume_index (vref.dim, i, j, k_p);
                    idx_n = volume_index (vref.dim, i, j, k_n);
                    out_img[gk] = (float)((ref_img[idx_n] - ref_img[idx_p]) / 2.0 / vref.pix_spacing[2]);
                }
            }
        }
        vout.img = out_img;
    }
    
    public int volume_index (int[] dims, int i, int j, int k)
    {
        return i + (dims[0] * (j + dims[1] * k));
    }

    public Volume vf_warp (Volume vout, Volume vin, Volume vf)
    {
        int d, i, j, k, v;
        int mi, mj, mk, mv;
        float fx, fy, fz;
        float mx, my, mz;
        float[] vf_img =  vf.img;
        float[] vin_img = vin.img;
        float[] vout_img;
        float[] dxyz = new float[3];
        if (vout == null) {
            vout = volume_create (vin.dim, vin.offset, vin.pix_spacing, vin.pix_type, vin.direction_cosines, 0);
        }
        vout_img = vout.img;

        /* Assumes size, spacing of vout same as size, spacing of vf */
        for (d = 0; d < 3; d++) {
            if (vout.dim[d] != vf.dim[d]) {
                System.err.printf("Dimension mismatch between fixed and moving\n");
                System.exit(1);
            }
            if (vout.pix_spacing[d] != vf.pix_spacing[d]) {
                System.err.printf("Resolutions mismatch between fixed and moving\n");
                System.exit(1);
            }
            if (vout.offset[d] != vf.offset[d]) {
                System.err.printf("offset mismatch between fixed and moving\n");
                System.exit(1);
            }
        }

        for (v = 0, k = 0, fz = vf.offset[2]; k < vf.dim[2]; k++, fz += vf.pix_spacing[2]) {
            for (j = 0, fy = vf.offset[1]; j < vf.dim[1]; j++, fy += vf.pix_spacing[1]) {
                for (i = 0, fx = vf.offset[0]; i < vf.dim[0]; i++, fx += vf.pix_spacing[0], v++) {
                    dxyz[0] = vf_img[3*v];
                    dxyz[1] = vf_img[3*v+1];
                    dxyz[2] = vf_img[3*v+2];
                    mz = fz + dxyz[2];
                    mk = round_int((mz - vin.offset[2]) / vin.pix_spacing[2]);

                    my = fy + dxyz[1];
                    mj = round_int((my - vin.offset[1]) / vin.pix_spacing[1]);

                    mx = fx + dxyz[0];
                    mi = round_int((mx - vin.offset[0]) / vin.pix_spacing[0]);

                    mv = (mk * vin.dim[1] + mj) * vin.dim[0] + mi;

                    if (mk < 0 || mk >= vin.dim[2]) continue;
                    if (mj < 0 || mj >= vin.dim[1]) continue;
                    if (mi < 0 || mi >= vin.dim[0]) continue;

                    vout_img[v] = vin_img[(int)mv];
                }
            }
        }
        vout.img = vout_img;
        return vout;
    }
    
    public int round_int(double x)
    {
        int rd;
        if((x)>=0)
        {
            rd = (int)((x)+0.5);
        }
        else
        {
            rd = (int)(-(-(x)+0.5));
        }
    return rd;
    }   
}
