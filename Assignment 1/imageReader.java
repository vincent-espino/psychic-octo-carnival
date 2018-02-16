
/**
 * Modified Java starter code for CSCI 576 - Homework 1, Programming Part.
 * Written by Vincent-Patrick Espino, USC ID 5002-2761-17.
 *
 * Assumptions made: all images will be of size 352 x 288.
 * According to Chloe Legendre, this is a reasonable assumption to make, as long as it is noted.
 */

import java.awt.*;
import java.awt.image.*;
import java.io.*;
import javax.swing.*;

public class imageReader {

    JFrame frame;
    JLabel lbIm1;
    JLabel lbIm2;
    BufferedImage img;
    BufferedImage modifiedImage;

    // According to the professor, we are allowed to assume these dimensions for all images.
    public static final int WIDTH = 352;
    public static final int HEIGHT = 288;

    // Matrix used for converting RGB values to their respective YUV values.
    public static final double[][] rgb_to_yuv_matrix = {
            {0.299, 0.587, 0.114}, // Y channel
            {0.596, -0.274, -0.322}, // U channel
            {0.211, -0.523, 0.312} // V channel
    };

    public static final double[][] yuv_to_rgb_matrix = {
            {1.000, 0.956, 0.621}, // R channel
            {1.000, -0.272, -0.647}, // G channel
            {1.000, -1.106, 1.703} // B channel
    };


    /**
     * Converts the given RGB triple in array B into its YUV representation, based on the provided
     * rgb_to_yuv_matrix. Assume this method works until other issues arise.
     * @param B, the RGB triple to be converted.
     * @return computed_matrix, the YUV triple.
     */
    public double[] convertToYUV(double[] B) {
        double[] computed_matrix = new double[3];
        int current_index = 0;
        for(int i = 0; i < rgb_to_yuv_matrix.length; i++) {
            for(int j = 0; j < rgb_to_yuv_matrix[i].length; j++) {
                computed_matrix[current_index] += rgb_to_yuv_matrix[i][j] * B[j];
            }
            current_index++;
        }
        return computed_matrix;
    }

    /**
     * Converts the given YUV triple in array B into its RGB representation, based on the provided
     * yuv_to_rgb_matrix.
     * @param B, the YUV triple to be converted.
     * @return computed_matrix, the RGB triple.
     */
    public double[] convertToRGB(double[] B) {
        double[] computed_matrix = new double[3];
        int current_index = 0;
        for(int i = 0; i < yuv_to_rgb_matrix.length; i++) {
            for(int j = 0; j < yuv_to_rgb_matrix[i].length; j++) {
                computed_matrix[current_index] += yuv_to_rgb_matrix[i][j] * B[j];
            }
            current_index++;
        }
        return computed_matrix;
    }


    /**
     * Subsamples the given yuv_matrix by the rate sampling_rate in the Y channel (luminance).
     * A similar method exists for the U and V channels. Note that this also performs the upsampling.
     * @param yuv_matrix, the YUV representation matrix of the image.
     * @param sampling_rate, the rate at which the specific channel should be sampled. Note that a sampling
     *                       rate means that, given a rate N, throw away every N'th signal.
     *                       In this method, this means that at every N'th entry, remove the Y channel (set it to 0?).
     *                       We may assume subsampling occurs on the width dimension (scan horizontally?).
     * @return yuv_matrix, the same matrix passed in, but subsampled.
     */
    public double[][][] subsampleInY(double[][][] yuv_matrix, int sampling_rate) {
        // If the sampling rate is 1, this implies no subsampling. Return the matrix as given.
        if(sampling_rate <= 1) {
            return yuv_matrix;
        }
        else {
            for(int y = 0; y < HEIGHT; y++) {
                // If our current index is an index we should subsample, do stuff!
                int previous_neighbor = 0;
                int next_neighbor = 0;
                for(int x = 0; x < WIDTH; x++) {
                    if(x % sampling_rate != 0) {
                        double[] yuv_temp = yuv_matrix[y][x];
                        yuv_temp[0] = (yuv_matrix[y][previous_neighbor][0] + yuv_matrix[y][next_neighbor][0])/2;
                    }
                    else {
                        // Get the closest available neighbor before and after the current pixel,
                        // using sampling rate to find the next neighbor.
                        previous_neighbor = x;
                        next_neighbor = previous_neighbor + sampling_rate;
                        if(next_neighbor >= WIDTH) {
                            next_neighbor = previous_neighbor;
                        }
                    }

                }
            }
            return yuv_matrix;
        }
    }

    public double[][][] subsampleInU(double[][][] yuv_matrix, int sampling_rate) {
        // If the sampling rate is 1, this implies no subsampling. Return the matrix as given.
        if(sampling_rate <= 1) {
            return yuv_matrix;
        }
        else {
            for(int y = 0; y < HEIGHT; y++) {
                // If our current index is an index we should subsample, do stuff!
                int previous_neighbor = 0;
                int next_neighbor = 0;

                for(int x = 0; x < WIDTH; x++) {
                    double[] yuv_temp = yuv_matrix[y][x];
                    if(x % sampling_rate != 0) {
                        yuv_temp[1] = (yuv_matrix[y][previous_neighbor][1] + yuv_matrix[y][next_neighbor][1])/2;
                        yuv_matrix[y][x] = yuv_temp;
                    }
                    else {
                        // Get the closest available neighbor before and after the current pixel,
                        // using sampling rate to find the next neighbor.
                        previous_neighbor = x;
                        next_neighbor = previous_neighbor + sampling_rate;
                        if(next_neighbor >= WIDTH) {
                            next_neighbor = previous_neighbor;
                        }
                    }

                }
            }
            return yuv_matrix;
        }
    }

    public double[][][] subsampleInV(double[][][] yuv_matrix, int sampling_rate) {
        // If the sampling rate is 1, this implies no subsampling. Return the matrix as given.
        if(sampling_rate <= 1) {
            return yuv_matrix;
        }
        else {
            for(int y = 0; y < HEIGHT; y++) {
                // If our current index is an index we should subsample, do stuff!
                int previous_neighbor = 0;
                int next_neighbor = 0;
                for(int x = 0; x < WIDTH; x++) {
                    if(x % sampling_rate != 0) {
                        double[] yuv_temp = yuv_matrix[y][x];
                        yuv_temp[2] = (yuv_matrix[y][previous_neighbor][2] + yuv_matrix[y][next_neighbor][2])/2;
                    }
                    else {
                        // Get the closest available neighbor before and after the current pixel,
                        // using sampling rate to find the next neighbor.
                        previous_neighbor = x;
                        next_neighbor = previous_neighbor + sampling_rate;
                        if(next_neighbor >= WIDTH) {
                            next_neighbor = previous_neighbor;
                        }
                    }

                }
            }
            return yuv_matrix;
        }
    }

    /**
     * Quantizes each RGB triple within rgb_matrix, using quantizationValue to determine
     * how each value is quantized. If the value is <= 0 or >= 256, the original matrix
     * is returned and no changes to it are made.
     * @param rgb_matrix, the RGB values to be quantized.
     * @param quantizationValue, the quantization value. This value should be from 0-255 (256 implies no subsampling).
     * @return
     */
    public double[][][] quantize(double[][][] rgb_matrix, double quantizationValue) {
        if(quantizationValue >= 256 || quantizationValue < 0) {
            System.out.println("No quantization needed.");
            return rgb_matrix;
        }
        else if(quantizationValue >= 0 && quantizationValue <= 255) {
            // This is used to start the calculations of every interval level after the first.
            double interval = -1.0;

            // Can't quantize with a value of 0... basically means no color.
            if(quantizationValue == 0) {
                quantizationValue = 1;
            }

            // An interval size is 256 / Q. So, if Q = 64, for example, there will be
            // 64 intervals, each of size 4 between them (so 0, 3, 7, etc.)
            double intervalSize = 256.0 / quantizationValue;
            double[] quantizationIntervals = new double[(int)quantizationValue];
            interval += intervalSize;

            // The first interval should always be 0.
            quantizationIntervals[0] = 0;
            for(int i = 1; i < quantizationIntervals.length; i++) {
                quantizationIntervals[i] = interval;
                interval += intervalSize;
            }

            for(int y = 0; y < HEIGHT; y++) {
                for(int x = 0; x < WIDTH; x++) {
                    double[] rgb_temp = rgb_matrix[y][x];
                    double[] quantizedRGB = new double[3];
                    for(int j = 0; j < rgb_temp.length; j++) {
                        double quantizedLevel = (rgb_temp[j] - (rgb_temp[j] % intervalSize))/intervalSize;
                        double quantizeToCheck = 0.0;
                        if(quantizedLevel >= quantizationIntervals.length - 1) {
                            quantizedRGB[j] = quantizationIntervals[(int)(quantizationIntervals.length - 1)];
                        }
                        else if(quantizedLevel <= 0) {
                            quantizedRGB[j] = quantizationIntervals[0];
                        }
                        else {
                            quantizeToCheck = (quantizationIntervals[(int)quantizedLevel] + quantizationIntervals[(int)quantizedLevel + 1])/2;
                            double quantizedValue = 0.0;
                            if(rgb_temp[j] < quantizeToCheck) {
                                quantizedValue = quantizationIntervals[(int)quantizedLevel];
                            }
                            else {
                                quantizedValue = quantizationIntervals[(int)quantizedLevel + 1];
                            }
                            quantizedRGB[j] = quantizedValue;
                        }
                    }
                    rgb_matrix[y][x] = quantizedRGB;
                }
            }
            return rgb_matrix;
        }
        return null;
    }

    public void ComputeError(double[][][] originalImageRepresentation, double[][][] modifiedImageRepresentation) {

        // Error is defined as the Euclidean distance between each r, g, b value.
        // The Euclidean distance equation used is inspired by the following Wikipedia article:
        // https://en.wikipedia.org/wiki/Color_quantization#Algorithms

        double totalError = 0.0;
        for(int y = 0; y < HEIGHT; y++) {
            for(int x = 0; x < WIDTH; x++) {
                double[] originalRGB = originalImageRepresentation[y][x];
                double[] modifiedRGB = modifiedImageRepresentation[y][x];
                double rDifference = (originalRGB[0] - modifiedRGB[0]) * (originalRGB[0] - modifiedRGB[0]);
                double gDifference = (originalRGB[1] - modifiedRGB[1]) * (originalRGB[1] - modifiedRGB[1]);
                double bDifference = (originalRGB[2] - modifiedRGB[2]) * (originalRGB[2] - modifiedRGB[2]);

                totalError += Math.sqrt(rDifference + gDifference + bDifference);
            }
        }

        System.out.println("Total error before dividing: " + totalError);
        totalError = totalError / (352.0 * 288.0);
        System.out.println("Error is " + totalError);
    }

    public void showIms(String[] args){
        if(args.length < 5) {
            System.out.println("Not enough command line arguments specified.");
            return;
        }

        // Parses the input from the command line, getting the Y, U, and V subsampling rates, as well as
        // the quantization value Q.
        int subsampling_Y = Integer.parseInt(args[1]);
        int subsampling_U = Integer.parseInt(args[2]);
        int subsampling_V = Integer.parseInt(args[3]);
        int quantization = Integer.parseInt(args[4]);
        System.out.format("Y: %d\nU: %d\nV: %d\nQ: %d\n", subsampling_Y, subsampling_U, subsampling_V, quantization);

        img = new BufferedImage(WIDTH, HEIGHT, BufferedImage.TYPE_INT_RGB);
        modifiedImage = new BufferedImage(WIDTH, HEIGHT, BufferedImage.TYPE_INT_RGB);
        double[][][] yuv_representation = new double[HEIGHT][WIDTH][];
        double[][][] original_rgb = new double[HEIGHT][WIDTH][];

        try {
            File file = new File(args[0]);
            InputStream is = new FileInputStream(file);

            long len = file.length();
            byte[] bytes = new byte[(int)len];

            int offset = 0;
            int numRead = 0;
            while (offset < bytes.length && (numRead=is.read(bytes, offset, bytes.length-offset)) >= 0) {
                offset += numRead;
            }


            int ind = 0;
            for(int y = 0; y < HEIGHT; y++){

                for(int x = 0; x < WIDTH; x++){

                    byte a = 0;
                    byte r = bytes[ind];
                    byte g = bytes[ind+HEIGHT*WIDTH];
                    byte b = bytes[ind+HEIGHT*WIDTH*2];

                    int pix = 0xff000000 | ((r & 0xff) << 16) | ((g & 0xff) << 8) | (b & 0xff);
                    double[] rgb_values = {(r & 0xff), (g & 0xff), (b & 0xff)};
                    original_rgb[y][x] = rgb_values;

                    // 1. Convert to YUV space.
                    double[] yuv_values = convertToYUV(rgb_values);
                    yuv_representation[y][x] = yuv_values;

                    //int pix = ((a << 24) + (r << 16) + (g << 8) + b);
                    img.setRGB(x,y,pix);

                    ind++;
                }
            }

            // 2. and 3. Process YUV subsampling and upsample for display.
            yuv_representation = subsampleInY(yuv_representation, subsampling_Y);
            yuv_representation = subsampleInU(yuv_representation, subsampling_U);
            yuv_representation = subsampleInV(yuv_representation, subsampling_V);

            // 4. Convert back to the RGB color space.
            for(int y = 0; y < HEIGHT; y++) {
                for(int x = 0; x < WIDTH; x++) {
                    double[] rgb_temp = yuv_representation[y][x];
                    rgb_temp = convertToRGB(rgb_temp);
                    yuv_representation[y][x] = rgb_temp;
                }
            }

            // 5. Quantize RGB channels.
            yuv_representation = quantize(yuv_representation, (double) quantization);


            for(int y = 0; y < HEIGHT; y++) {
                for(int x = 0; x < WIDTH; x++) {
                    double[] yuv_temp = yuv_representation[y][x];

                    // Fixes negative/too large of values, since RGB can only range from 0-255.
                    if(yuv_temp[0] < 0) {
                        yuv_temp[0] = 0;
                    }
                    else if(yuv_temp[0] > 255) {
                        yuv_temp[0] = 255;
                    }
                    if(yuv_temp[1] < 0) {
                        yuv_temp[1] = 0;
                    }
                    else if(yuv_temp[1] > 255) {
                        yuv_temp[1] = 255;
                    }
                    if(yuv_temp[2] < 0) {
                        yuv_temp[2] = 0;
                    }
                    else if(yuv_temp[2] > 255) {
                        yuv_temp[2] = 255;
                    }
                    int yuv_pix = 0xff000000 | (((int)yuv_temp[0] & 0xff) << 16) | (((int)yuv_temp[1] & 0xff) << 8) | ((int)yuv_temp[2] & 0xff);
                    modifiedImage.setRGB(x, y, yuv_pix);
                }
            }


//            System.out.println("Processing done.");

//            ComputeError(original_rgb, yuv_representation);
//            System.out.format("Error computed with values of Y: %d, U: %d, V: %d, Q: %d\n", subsampling_Y, subsampling_U, subsampling_V, quantization);

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }

        // Use labels to display the images
        frame = new JFrame();
        GridBagLayout gLayout = new GridBagLayout();
        frame.getContentPane().setLayout(gLayout);

        JLabel lbText1 = new JLabel("Original image (Left)");
        lbText1.setHorizontalAlignment(SwingConstants.CENTER);
        JLabel lbText2 = new JLabel("Image after modification (Right)");
        lbText2.setHorizontalAlignment(SwingConstants.CENTER);
        lbIm1 = new JLabel(new ImageIcon(img));
        lbIm2 = new JLabel(new ImageIcon(modifiedImage));

        GridBagConstraints c = new GridBagConstraints();
        c.fill = GridBagConstraints.HORIZONTAL;
        c.anchor = GridBagConstraints.CENTER;
        c.weightx = 0.5;
        c.gridx = 0;
        c.gridy = 0;
        frame.getContentPane().add(lbText1, c);

        c.fill = GridBagConstraints.HORIZONTAL;
        c.anchor = GridBagConstraints.CENTER;
        c.weightx = 0.5;
        c.gridx = 1;
        c.gridy = 0;
        frame.getContentPane().add(lbText2, c);

        c.fill = GridBagConstraints.HORIZONTAL;
        c.gridx = 0;
        c.gridy = 1;
        frame.getContentPane().add(lbIm1, c);

        c.fill = GridBagConstraints.HORIZONTAL;
        c.gridx = 1;
        c.gridy = 1;
        frame.getContentPane().add(lbIm2, c);

        frame.pack();
        frame.setVisible(true);

    }

    public static void main(String[] args) {
        imageReader ren = new imageReader();
        ren.showIms(args);
    }

}