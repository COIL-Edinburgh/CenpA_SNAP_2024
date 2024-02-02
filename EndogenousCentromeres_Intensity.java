import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.Roi;
import ij.io.DirectoryChooser;
import ij.plugin.ChannelSplitter;
import ij.plugin.ImageCalculator;
import ij.plugin.PlugIn;
import ij.plugin.frame.RoiManager;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;

import java.awt.*;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Paths;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;

public class EndogenousCentromeres_Intensity implements PlugIn {

    boolean changeParameters;
    String baseDirectory;
    String filename;
    double squareSize;
    double minCircularity;
    double maxFeret;
    double minCentromere;
    double maxCentromere;
    double caHor;
    double caVer;
    int slices;
    double scale;
    RoiManager roiManager;
    String filePath;

    public void run(String arg) {

        DirectoryChooser dc = new DirectoryChooser("Bio-Formats Mass Importer");
        baseDirectory = dc.getDirectory();

        // list of files to actually open with Bio-Formats Importer
        ArrayList<String> filesToOpen = new ArrayList<>();

        // process all files in the chosen directory
        File dir = new File(baseDirectory);
        File[] files = dir.listFiles();
        boolean firstFile=true;
        roiManager = new RoiManager();

        //create a new folder to save files in
        String newDirectory = baseDirectory + "Output_v2_4Channels";
        filePath = Paths.get(baseDirectory, "Output_v2_4Channels").toString();
        new File(newDirectory).mkdir();
        int[] fileNameList = new int[3]; //To hold the channel order

        //For each file
        for (int m = 0; m < files.length; m++) {
            filesToOpen.add(files[m].getPath());
            String id = filesToOpen.get(m);

            //Check it is an .ims of .tif
            if ((id.contains(".ims")||id.contains(".tif"))&& (!id.contains(".xml") && !id.contains(".zip"))) {

                //Open the file get info and split channels, make z-projections
                IJ.run("Bio-Formats Importer", "open=" + id + " color_mode=Default open_files view=Hyperstack stack_order=XYCZT");
                ImagePlus imp = WindowManager.getCurrentImage();
                filename = imp.getShortTitle();
                slices = WindowManager.getCurrentImage().getNSlices();
                scale = imp.getCalibration().pixelWidth;
                new ChannelSplitter();
                ImagePlus[] channels = ChannelSplitter.split(imp);
                if (channels.length == 4) {
                    for (int i = 0; i < channels.length; i++) {
                        IJ.run(channels[i], "Z Project...", "projection=[Max Intensity]");
                        WindowManager.getCurrentImage().show();
                        WindowManager.getCurrentImage().setTitle(i + "_Zproj");
                    }

                    //If it is the first image file set the channel order and parameters and make the logfile
                    if (firstFile) {
                        fileNameList = channelSelector(channels.length);
                        setParameterSettings();
                        makeLogFile(fileNameList);
                        firstFile = false;
                    }

                    //Make a composite file of the 4 channels
                    IJ.run("Merge Channels...", "c1=[" + fileNameList[0] + "_Zproj] c2=[" + fileNameList[1] + "_Zproj] c3=[" + fileNameList[2] + "_Zproj] c7=[" + fileNameList[3] + "_Zproj] create keep");
                    ImagePlus compImage = WindowManager.getCurrentImage();

                    //Make the DAPI and Reference image masks
                    ImagePlus maskDAPI = makeDAPImask(WindowManager.getImage(fileNameList[2] + "_Zproj"));
                    //ImagePlus maskRef = makeRefMask(WindowManager.getImage(fileNameList[1] + "_Zproj"), maskDAPI);
                    ImagePlus maskRef = WindowManager.getImage(fileNameList[1] + "_Zproj");

                    //Get Outlines of nucleii, label them and draw onto the composite image, save the image
                    Roi[] cellRois = getCellRois(maskDAPI);
                    roiManager.reset();

                    Roi[] rois = getCentromereRois(maskRef, cellRois);

                    //From the reference image mask find the ROIs that fit the criteria from the parameter settings
                    maskRef.setTitle("maskref");
                    maskRef.show();
                    IJ.saveAs("Tiff", Paths.get(filePath, filename + "_mask_v2_4Channels.tif").toString());

                    //Draw boxes on the data z project and get the intensities
                    double[][] data = getData(rois, WindowManager.getImage(fileNameList[0] + "_Zproj"), WindowManager.getImage(fileNameList[3]+"_Zproj"), cellRois);

                    //Save the ROIs and output the data to the data file
                    roiManager.runCommand("Save", Paths.get(filePath, filename + "_v2_4Channels.zip").toString());
                    outputData(data);

                    compImage = drawCells(compImage, cellRois, rois);
                    IJ.saveAs("Tiff", Paths.get(filePath, filename + "_v2_4Channels.tif").toString());
                    IJ.log(Paths.get(filePath, filename + "_v2_4Channels.tif").toString());
                }else {
                    IJ.log(imp.getShortTitle()+ " does not have 4 channels and will not be processed");
                }
                //Reset roi manager and close all windows
                roiManager.reset();
                IJ.run("Close All", "");


            }
        }
        roiManager.close();
    }


    private Roi[] getCentromereRois(ImagePlus mask, Roi[] cells){
        for(Roi cell:cells){
            mask.setRoi(cell);
            IJ.setAutoThreshold(mask, "Triangle dark");
            IJ.run(mask, "Analyze Particles...", "size=" + minCentromere + "-" + maxCentromere + " pixel circularity=" + minCircularity + "-1.00 show=Nothing add");
        }
        return roiManager.getRoisAsArray();
    }

    //Opens a dialog box to set the Reference, Data and DAPi channels and select/deselect default parameter settings.
    // Takes an interger number of channels in the input image and returns an integer array for the channel order.
    public int[] channelSelector(int channels){
        GenericDialog channelDialog = new NonBlockingGenericDialog("Channel Selector");
        int[] filenameList = new int[channels];
        for (int k = 0; k< channels; k++){filenameList[k]= k;}
        String[] coloursArray = new String[]{"Data Channel","Reference Channel","DAPI Channel","Other Channel?"};
        channelDialog.addNumericField(coloursArray[0], 2);
        channelDialog.addNumericField(coloursArray[1], 3);
        channelDialog.addNumericField(coloursArray[2], 0);
        channelDialog.addNumericField(coloursArray[3], 1);
        channelDialog.addCheckbox("Change default parameter settings?", false );
        channelDialog.showDialog();
        int[] choicesArray = new int[4];
        for (int i = 0; i < 4; i++){
            choicesArray[i] = (int) channelDialog.getNextNumber();
        }
        changeParameters= channelDialog.getNextBoolean();

        if ((choicesArray[0]==choicesArray[2]||choicesArray[1]==choicesArray[2])){
            System.out.print("Reference and Data channels should be different from DAPI channel");
        }else if(choicesArray[0]==choicesArray[1]) {
            IJ.log("WARNING! Reference and Data channels are the same");
        }
        return choicesArray;
    }

    private void makeLogFile(int[] channels){

        Date date = new Date(); // This object contains the current date value
        SimpleDateFormat formatter = new SimpleDateFormat("dd-MM-yyyy HH:mm:ss");
        String CreateName = Paths.get( filePath, "Logfile_OUTPUT.txt").toString();
        IJ.log(CreateName);
        try{
            FileWriter fileWriter = new FileWriter(CreateName,true);
            BufferedWriter bufferedWriter = new BufferedWriter(fileWriter);
            bufferedWriter.newLine();
            bufferedWriter.write(formatter.format(date));
            bufferedWriter.newLine();
            bufferedWriter.write("Method: Bandpass; minmax");
            bufferedWriter.newLine();
            bufferedWriter.write("CRaQ_Toni_Plugin version");
            bufferedWriter.newLine();
            bufferedWriter.write("Base Directory: "+ baseDirectory);
            bufferedWriter.newLine();
            bufferedWriter.write("Reference Channel: "+channels[1]);
            bufferedWriter.newLine();
            bufferedWriter.write("Data Channel: "+channels[0]);
            bufferedWriter.newLine();
            bufferedWriter.write("DAPI Channel: "+channels[2]);
            bufferedWriter.newLine();
            bufferedWriter.write("Other Channel: "+channels[3]);
            bufferedWriter.newLine();
            bufferedWriter.newLine();
            bufferedWriter.write ("Square Size: "+ squareSize);
            bufferedWriter.newLine();
            bufferedWriter.write("Minimum Circularity: "+minCircularity);
            bufferedWriter.newLine();
            bufferedWriter.write("Maximum Ferets Diameter: "+ maxFeret);
            bufferedWriter.newLine();
            bufferedWriter.write ("Minimum Centromere Size: "+ minCentromere);
            bufferedWriter.newLine();
            bufferedWriter.write ("Maximum Centromere Size: "+maxCentromere);
            bufferedWriter.newLine();
            bufferedWriter.write("Chromatic aberration correction: ("+caHor+","+caVer+") [(x,y) difference of reference compared to data]");
            bufferedWriter.newLine();
            bufferedWriter.newLine();
            bufferedWriter.close();
        }
        catch(IOException ex) {
            System.out.println(
                    "Error writing to file '" + CreateName + "'");
        }
    }

    //Opens the parameter settings dialog box and allows users to input values for variable parameter settings, updates
    // the global variables
    public void setParameterSettings(){
        GenericDialog parameterDialog = new NonBlockingGenericDialog("Change parameter Settings");
        parameterDialog.addNumericField("Square size",7,0,0,"pixels");
        parameterDialog.addNumericField("Minimum Circularity",0.95,2,4,"a.u.");
        parameterDialog.addNumericField("Max Feret's Diameter",7,1,3,"pixels");
        parameterDialog.addNumericField("Min Centromere Size",4,0,2,"pixel");
        parameterDialog.addNumericField("Max Centromere Size",35,0,2,"pixel");
        parameterDialog.addMessage("\nIf known, set the chromatic aberration of the reference channel compared to the data channel.");
        parameterDialog.addNumericField("Chromatic aberration (horizontal): ",0,0,2,"pixels to right");
        parameterDialog.addNumericField("Chromatic aberration (vertical): ",0,0,2,"pixels down");
        if (changeParameters) {
            parameterDialog.showDialog();
        }
        squareSize = parameterDialog.getNextNumber();
        minCircularity = parameterDialog.getNextNumber();
        maxFeret = parameterDialog.getNextNumber();
        minCentromere = parameterDialog.getNextNumber();
        maxCentromere = parameterDialog.getNextNumber();

        caHor = parameterDialog.getNextNumber();
        caVer = parameterDialog.getNextNumber();

        if (minCircularity>=1){
            System.out.print("Minimum Circularity should be smaller than 1");
        }
        if (minCentromere>=maxCentromere){
            System.out.print("Minimum Centromere size should be smaller than maximum centromere size");
        }

    }

    //Takes the ImagePlus for the Z-projected dapi channel and returns a nucleii mask
    private ImagePlus makeDAPImask(ImagePlus dapi) {
        WindowManager.getImage(dapi.getShortTitle());
        IJ.run(dapi,"Duplicate...", "title=blur");
        IJ.selectWindow("blur");
        ImagePlus blur = WindowManager.getCurrentImage();
        IJ.run("Gaussian Blur...", "sigma=75");

        IJ.run(dapi, "Duplicate...", "title=blur5");
        IJ.selectWindow("blur5");
        ImagePlus blur5 = WindowManager.getCurrentImage();
        IJ.run("Gaussian Blur...", "sigma=5");

        IJ.setAutoThreshold(blur5, "Huang dark");
        IJ.run(blur5, "Analyze Particles...", "show=Masks");
        IJ.run("Watershed");
        IJ.run("Erode");
        IJ.run( "Analyze Particles...", "size= 50-Infinity show=Masks");
        //IJ.run("Invert");
        IJ.run("16-bit");
        IJ.run("Multiply...", "value=257.000");
        IJ.run("Invert");
        ImagePlus mask = WindowManager.getCurrentImage();
        blur.changes=false;
        blur.close();
        blur5.changes=false;
        blur5.close();
        return mask;
    }


    //Takes the Dapi mask image and applies a threshold then uses analyse particles to return ROI outlines of the nucleii
    public Roi[] getCellRois(ImagePlus Dapi){
        int maskDAPIID = Dapi.getID();
        IJ.selectWindow(maskDAPIID);
        IJ.setAutoThreshold(Dapi, "Default dark");//To work to my computer use Default dark???
        IJ.selectWindow(maskDAPIID);
        IJ.run(Dapi, "Analyze Particles...", "size=5000-Infinity pixel add");
        return roiManager.getRoisAsArray();
    }

    //Draws and labels the Nucleii outlines in the Roi[] input onto the comp image input and returns a new composite image.
    private ImagePlus drawCells(ImagePlus comp, Roi[] cells, Roi[] centromeres){
        ImagePlus[] channels = ChannelSplitter.split(comp);
        for (int j=0;j<4;j++){
            channels[j].show();}
        ImagePlus draw = IJ.createImage("draw", channels[0].getWidth(), channels[0].getHeight(), 1, channels[0].getBitDepth());
        draw.show();
        for(int i = 0; i < cells.length; i++){
            IJ.setForegroundColor(255, 255, 255);
            ImageProcessor ip = draw.getProcessor();
            Font font = new Font("SansSerif", Font.BOLD, 20);
            ip.setFont(font);
            ip.setColor(Color.white);
            String cellnumber = String.valueOf(i+1);
            int xpos = (int) cells[i].getContourCentroid()[0];
            int ypos = (int) cells[i].getContourCentroid()[1];
            ip.drawString(cellnumber, xpos, ypos);
            ip.draw(cells[i]);
            draw.updateAndDraw();
        }
        for(int j = 0; j < centromeres.length; j++){
            IJ.setForegroundColor(255, 255, 255);
            ImageProcessor ip = draw.getProcessor();
            Font font = new Font("SansSerif", Font.BOLD, 10);
            ip.setFont(font);
            ip.setColor(Color.white);
            int xpos = (int) centromeres[j].getContourCentroid()[0];
            int ypos = (int) centromeres[j].getContourCentroid()[1];
            ip.draw(centromeres[j]);
            draw.updateAndDraw();
        }

        IJ.run( "Merge Channels...", "c1="+channels[0].getShortTitle()+" c2="+channels[1].getShortTitle()+"" +
                " c3="+channels[2].getShortTitle()+" c4="+draw.getShortTitle()+" c7="+channels[3].getShortTitle()+" create keep");
        return WindowManager.getCurrentImage();
    }

    //Takes the Roi[] rois which specifies the thresholded spots in the reference channel and applies them to the data
    // channel. Returns a double[3][roi.length] with the roi number, intensity in the data channel and which nucleii (cellRois)
    //it is contain in for each non-overlapping roi.
    private double[][] getData(Roi[] rois, ImagePlus dataZproj, ImagePlus otherZproj, Roi[] cellRois) {
        roiManager.reset();
        double corner = (squareSize - 1) / 2;
        int count = 1;
        double[][] data= new double[4][rois.length];
        for (int i = 0; i < rois.length; i++) {
            if (rois[i].getFeretsDiameter() < maxFeret) {
                double x = rois[i].getContourCentroid()[0];
                double y = rois[i].getContourCentroid()[1];
                double cx = x - caHor;
                double cy = y - caVer;
                Roi rectangle = new Roi(cx - corner, cy - corner, squareSize, squareSize);
                dataZproj.setRoi(rectangle);
                ImageStatistics stats = dataZproj.getStatistics();
                otherZproj.setRoi(rectangle);
                ImageStatistics otherStats = otherZproj.getStatistics();
                if (stats.min > 0 && stats.max < 65000 && stats.area == (squareSize * squareSize * scale * scale)) {
                    if (stats.max > 0) {
                        data[0][i]=count;
                        data[1][i]=stats.max-stats.min;
                        data[2][i]=whichCell(x,y,cellRois);
                        data[3][i]=otherStats.max-otherStats.min;
                        roiManager.addRoi(rectangle);
                    } else {
                        data[0][i]=count;
                    }
                    IJ.run(dataZproj, "Clear", "slice");
                    count++;
                }
            }
        }
        return data;
    }

    private int whichCell(double x, double y, Roi[] cellRois){
        int cell = 0;
        for(int i = 0; i < cellRois.length; i++){
            if(cellRois[i].contains((int)x,(int)y)){
                cell=i+1;
            }
        }
        return cell;
    }

    private void outputData(double[][] data){

        Date date = new Date(); // This object contains the current date value
        SimpleDateFormat formatter = new SimpleDateFormat("dd-MM-yyyy , HH:mm:ss");
        String CreateName = Paths.get(filePath, "OUTPUT.csv").toString();

        try{
            FileWriter fileWriter = new FileWriter(CreateName,true);
            BufferedWriter bufferedWriter = new BufferedWriter(fileWriter);
            bufferedWriter.newLine();
            bufferedWriter.write(formatter.format(date));
            bufferedWriter.newLine();
            bufferedWriter.write(filename);
            bufferedWriter.newLine();
            bufferedWriter.write("N, Max-min Data, Max-min Other, Cell");
            bufferedWriter.newLine();
            for (int i =0; i< data[0].length; i++){
                if(data[0][i]!=0 && data[2][i]!=0){
                    bufferedWriter.write(data[0][i]+","+data[1][i]+","+data[3][i]+ ","+ data[2][i]);
                    bufferedWriter.newLine();
                }}
            bufferedWriter.close();
        }
        catch(IOException ex) {
            System.out.println(
                    "Error writing to file '" + CreateName + "'");
        }
    }


}