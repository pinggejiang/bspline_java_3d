package Bspline_Java_3D;
import java.util.Arrays;

/*
 *
 * @author pingge
 * May 20th, 2014
 * File usage: class definition
 *
 */

public class Volume {
	    public int[] dim = new int[3];			    // x, y, z Dims
	    public int npix;			                // # of voxels in volume
	    public float[] offset = new float[3];
	    public float[] pix_spacing = new float[3];	// voxel spacing
	    public float[] direction_cosines = new float[9];
        public String pix_type;	                    // Voxel Data type
	    public int pix_size;		                // # bytes per voxel
	    public float[] img;			                // Voxel Data
	
}
