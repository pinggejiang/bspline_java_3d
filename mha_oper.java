package Bspline_Java_3D;
import java.io.*;
import java.util.*;
import java.nio.channels.FileChannel;
import java.nio.MappedByteBuffer;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
/*
 *
 * @author pingge
 * May 20th, 2014
 * File usage: function definition, for read_mha and write_mha
 *
 */

public class mha_oper {
    /*different image file extensions, check if mh5*/
    public static boolean is_mh5 (String filename)
    {
        int len = filename.length();
        if (len < 4 )
        {
           return false;
        }
        return filename.substring(len-4).equals(".mh5") || filename.substring(len-4).equals(".MH5");
        
    }
    /*read image files based on extension*/
    public Volume Read_mha (String filename)
    {
        if(is_mh5(filename))
        {
            return read_mha_internal(filename, 1);
        }
        else
        {
            return read_mha_internal(filename, 0);
        }
    }
    /*read image files main*/
    public Volume read_mha_internal(String filename, int mh5)
    {
        try
        {
            Volume vol = new Volume();
            System.out.printf("Reading file: %s\n", filename);
            FileReader fr = new FileReader (filename);
            BufferedReader br = new BufferedReader(fr);
            String linebuff;
            vol.pix_type = "PT_UNDEFINED";
            vol.pix_size = -1;
            while ((linebuff = br.readLine()) !=null)
            {
                String[] token = linebuff.split(" ");
                if (token[1].equals("="))
                {
                    if (token[0].equals("ElementDataFile"))
                    {
                        if (!token[2].equals("LOCAL"))
                        {
                            break;
                        }
                        continue;
                    }

                    if (token[0].equals("DimSize"))
                    {
                        vol.dim[0] = Integer.parseInt(token[2]);
                        vol.dim[1] = Integer.parseInt(token[3]);
                        vol.dim[2] = Integer.parseInt(token[4]);
                        vol.npix = vol.dim[0] * vol.dim[1] * vol.dim[2];
                        continue;
                    }
                    if (token[0].equals("Offset"))
                    {
                        vol.offset[0] = Float.parseFloat(token[2]);
                        vol.offset[1] = Float.parseFloat(token[3]);
                        vol.offset[2] = Float.parseFloat(token[4]);
                        continue;
                    }
                    if (token[0].equals("ElementSpacing"))
                    {
                        vol.pix_spacing[0] = Float.parseFloat(token[2]);
                        vol.pix_spacing[1] = Float.parseFloat(token[3]);
                        vol.pix_spacing[2] = Float.parseFloat(token[4]);
                        continue;
                    }
                    if (token[0].equals("ElementNumberOfChannels"))
                    {
                        if (vol.pix_type.equals("PT_UNDEFINED") || vol.pix_type.equals("PT_FLOAT"))
                        {
                            vol.pix_type = "PT_VF_FLOAT_INTERLEAVED";
                            vol.pix_size = 3 * 4; //4 is the size of float 
                        }
                        continue;
                    }
                    if (linebuff.equals("ElementType = MET_FLOAT"))
                    {
                        if (vol.pix_type.equals("PT_UNDEFINED"))
                        {
                            vol.pix_type = "PT_FLOAT";
                            vol.pix_size = 4;
                        }
                        continue;
                    }
                    if (linebuff.equals("ElementType = MET_SHORT"))
                    {
                        vol.pix_type = "PT_SHORT";
                        vol.pix_size = 2;
                        continue;
                    }
                    if (linebuff.equals("ElementType = MET_UCHAR"))
                    {
                        vol.pix_type = "PT_UCHAR";
                        vol.pix_size = 1;
                        continue;
                    }
                }
                else
                {
                    break;
                }
            }
            vol.img = new float[vol.npix];
            freadjava(vol.img, vol.pix_size, vol.npix, vol.pix_type, filename);
            vol.pix_size = 4;
            vol.pix_type = "PT_FLOAT";
            fr.close();
            return vol;
        }
        catch(Exception e)
        {
            System.out.printf("Read file failed\n", filename); 
            e.getMessage();
            return null;
        }
    }
    @SuppressWarnings("resource")
	public void freadjava(float[] img, int size, int npix, String ptype, String filename)
    {
        FileChannel inChannel;
        try {
            inChannel = new RandomAccessFile(filename, "r").getChannel();
            MappedByteBuffer buffer = inChannel.map(FileChannel.MapMode.READ_ONLY, 0, inChannel.size());
            byte[] bytes = new byte[(int)inChannel.size()];
            buffer.get(bytes, 0, bytes.length);
            byte[] img_bytes = new byte[npix * size];
            img_bytes = Arrays.copyOfRange(bytes, (int)inChannel.size() - size * npix , (int)inChannel.size() );
            switch (ptype) 
            {
                case "PT_SHORT":
                    short[] img_short = new short[npix];
                    ByteBuffer.wrap(img_bytes).order(ByteOrder.LITTLE_ENDIAN).asShortBuffer().get(img_short);
                    for (int i = 0; i < npix; i++)
                    {
                    	img[i] = (float)img_short[i];
                    }
                    break;
                case "PT_FLOAT":
                    ByteBuffer.wrap(img_bytes).order(ByteOrder.LITTLE_ENDIAN).asFloatBuffer().get(img);
                    break;
                default:
                    System.err.println("Currently unsupported binary type");
                    System.exit(1);
            }
            inChannel.close();
        }
        catch (Exception e)
        {
            System.err.println("Error reading image");
            System.exit(1);
        }
    }
    /*read image files based on extension*/
    public void write_mha (String filename, Volume vol)
    {
        if (is_mh5 (filename)) {
            write_mha_internal (filename, vol, 1);
        } else {
            write_mha_internal (filename, vol, 0);
        }
    }
    
    /*write mha main*/
    public void write_mha_internal (String filename, Volume vol, int mh5)
    {
        FileOutputStream out;
        PrintStream p;
        String element_type = null;
        try
        {
            out = new FileOutputStream(filename);
            p = new PrintStream(out);
            p.println("ObjectType = Image");
            p.println("NDims = 3");
            p.println("BinaryData = True");
            p.println("BinaryDataByteOrderMSB = False");
            p.println("TransformMatrix = 1 0 0 0 1 0 0 0 1");
            p.printf("Offset = %g %g %g\n", vol.offset[0], vol.offset[1], vol.offset[2]);
            p.println("CenterOfRotation = 0 0 0");
            p.printf("ElementSpacing = %g %g %g\n", vol.pix_spacing[0], vol.pix_spacing[1], vol.pix_spacing[2]);
            p.printf("DimSize = %d %d %d\n", vol.dim[0], vol.dim[1], vol.dim[2]);
            p.println("AnatomicalOrientation = RAI");
            switch(vol.pix_type)
            {
                case "PT_SHORT":
                    element_type = "MET_SHORT";
                    p.printf("ElementType = %s\n", element_type);
                    p.println("ElementDataFile = LOCAL");
                    ByteBuffer buffer = ByteBuffer.allocate(vol.npix * vol.pix_size);
                    buffer.order(ByteOrder.LITTLE_ENDIAN);
                    short[] srt_pix = new short[vol.npix];
                    for (int i = 0; i < vol.npix; i++)
                    {
                    	srt_pix[i] = (short)vol.img[i];
                    }
                    buffer.asShortBuffer().put(srt_pix);
                    byte[] bytes = buffer.array();
                    p.write(bytes, 0, vol.npix * vol.pix_size);
                    break;
                case "PT_FLOAT":
                    element_type = "MET_FLOAT";
                    p.printf("ElementType = %s\n", element_type);
                    p.println("ElementDataFile = LOCAL");
                    out.flush();
                    out.close();
                    byte[] btes = new byte[vol.npix * vol.pix_size];
                    int k = 0;
                    for (int i = 0; i < vol.npix; i++)
                    {
                    	int data = Float.floatToIntBits(vol.img[i]);
                    	btes[k++]=(byte)(data & 0xff);
                        btes[k++]=(byte)((data >> 8) & 0xff);
                        btes[k++]=(byte)((data >> 16) & 0xff);
                        btes[k++]=(byte)((data >> 24) & 0xff);
                    }
                    File f = new File(filename);
                    long fileLength = f.length();
                    RandomAccessFile raf = new RandomAccessFile(f, "rw");
                    raf.seek(fileLength);
                    raf.write(btes);
                    break;
                    default:
                    	System.out.println("Type not supported yet");
                    	System.exit(-1);
            }
            System.out.println("Done");
        }
        catch (Exception e)
        {
            System.err.println("error writing to file");
            System.exit(-1);
        }
    }

}
    

