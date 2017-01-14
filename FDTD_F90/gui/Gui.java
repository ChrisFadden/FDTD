/** @file Gui.java
 *  @brief GUI implementation for electromagnetic simulations
 *
 *  @author Chris Fadden
 *
 */

import java.awt.*;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import javax.swing.*;

import java.io.File;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

public class Gui extends JFrame implements ActionListener{
  
  private JTextField NlambdaField, TFSF_angleField;
  private JTextField TFSF_x0field, TFSF_y0field;
  private JTextField SizeXfield, SizeYfield;
  private JTextField maxTimeField, FreqField;
  private JTextField TFSF_sizeXfield, TFSF_sizeYfield;
  private JTextField PMLfield;
 
  private JRadioButton Harmonic, Gaussian;
  private JRadioButton TFSF, SoftSource;
  private ButtonGroup Source;
  private ButtonGroup SourceType; 

  private GridBagConstraints c; 
  static public JPanel panel;  
  private JButton doneButton; 
  
  private String SizeX, SizeY;
  private String Nlambda, maxTime;
  private String freq, TFSFangle;
  private String srcType, src;
  private String TFSF_x0, TFSF_y0;
  private String TFSF_sizeX, TFSF_sizeY;
  private String PML_Size; 

  private int row;

  /** @brief Constructor the GUI class
   *
   *    This uses the GridBagLayout to arrange
   *  textboxes for input to electromagnetic simulations.
   *  Default values are set, and the user can choose to 
   *  modify whichever parameters are preferred. 
   */
 public Gui() {
    
    panel = new JPanel();     
    panel.setLayout(new GridBagLayout());
    c = new GridBagConstraints(); 
    c.weighty = 0.5; 
    row = 0;
    
    /****************
     * Grid
     ***************/

    c.gridx = 0;
    c.gridy = row;
    c.ipady = 20; 
    panel.add(new JLabel("Grid"), c);
    row++;
    c.ipady = 0;
    
    SizeXfield = new JTextField("200     ");
    placeField("Size X :    ", SizeXfield);
    
    SizeYfield = new JTextField("200    ");
    placeField("Size Y :    ", SizeYfield);
      
    maxTimeField = new JTextField("500         ");  
    placeField("Max Time : ", maxTimeField);

    FreqField = new JTextField("1        ");  
    placeField("Frequency (GHz) : ", FreqField); 
 
    NlambdaField = new JTextField("50        ");  
    placeField("Points Per Wavelength : ", NlambdaField); 
    
    TFSF_x0field = new JTextField("20        ");  
    placeField("TF/SF x start : ", TFSF_x0field); 
    
    TFSF_sizeXfield = new JTextField("150       ");  
    placeField("TF/SF Size (X) : ", TFSF_sizeXfield); 
    
    TFSF_y0field = new JTextField("20        ");  
    placeField("TF/SF y start : ", TFSF_y0field); 
    
    TFSF_sizeYfield = new JTextField("150       ");  
    placeField("TF/SF Size (Y) : ", TFSF_sizeYfield); 
    
    TFSF_angleField = new JTextField("45        ");  
    placeField("Source Angle of Incidence (degrees) : ", TFSF_angleField); 
    
    PMLfield = new JTextField("10       ");  
    placeField("Size of PML : ", PMLfield); 
       
    TFSF = new JRadioButton("Total Field/Scattered Field");
    SoftSource = new JRadioButton("Soft Source");
    SourceType = new ButtonGroup();
    placeTFSF_Field();
  
    Harmonic = new JRadioButton("Harmonic");
    Gaussian = new JRadioButton("Gaussian");
    Source = new ButtonGroup();
    placeSourceField();
  
    /*****************
     * Done Button
     ****************/
    c.gridx = 0;
    c.gridy = row;
    doneButton = new JButton("Done"); 
    panel.add(doneButton, c);
    row++;
    
    doneButton.addActionListener(this);

  }//end constructor 
  
  /** @brief Done button event handler
   *
   *   When the Done button is pushed,
   *  the current value in the textbox
   *  is taken, and then printed to a file.
   */
  public void actionPerformed(ActionEvent event)
  {
    Object source = event.getSource();
    if(source == doneButton)
    {
      
      SizeX = SizeXfield.getText();
      SizeY = SizeYfield.getText();
      
      maxTime = maxTimeField.getText();
      freq = FreqField.getText();
      Nlambda = NlambdaField.getText(); 

      TFSFangle = TFSF_angleField.getText();
      TFSF_x0 = TFSF_x0field.getText();
      TFSF_y0 = TFSF_y0field.getText();
      TFSF_sizeX = TFSF_sizeXfield.getText();
      TFSF_sizeY = TFSF_sizeYfield.getText();
      PML_Size = PMLfield.getText();

      srcType = "1"; //Soft Source
      if(TFSF.isSelected())
        srcType = "0";

      src = "1"; //Gaussian Source
      if(Harmonic.isSelected())
        src = "0";
      
      printToFile();
      
    }
  }//end button press
  
  public void placeField(String label, JTextField field){
  
    c.gridx = 0;
    c.gridy = row;
    panel.add(new JLabel(label), c);
  
    c.gridx = 1;
    c.gridy = row;
    panel.add(field, c);
    row++;
   
    return;
  }

  public void placeSourceField(){
  
    c.gridx = 0;
    c.gridy = row;
    c.ipady = 20; 
    panel.add(new JLabel("Source:"), c);
    row++;
    c.ipady = 0;
    
    c.gridx = 0;
    c.gridy = row; 
    Source.add(Harmonic);
    Source.add(Gaussian);
  
    Harmonic.setSelected(true);

    panel.add(Harmonic, c);
    row++;
    c.gridy = row;
    panel.add(Gaussian, c);
    row++;
  
    return; 
  }

 public void placeTFSF_Field(){
  
    c.gridx = 0;
    c.gridy = row;
    c.ipady = 20; 
    panel.add(new JLabel("Source Type:"), c);
    row++;
    c.ipady = 0;
    
    c.gridx = 0;
    c.gridy = row; 
    SourceType.add(TFSF);
    SourceType.add(SoftSource);
  
    TFSF.setSelected(true);

    panel.add(TFSF, c);
    row++;
    c.gridy = row;
    panel.add(SoftSource, c);
    row++;
  
    return; 
  }

  public void printToFile(){
     
    try{

      FileWriter writer = new FileWriter("grid.txt");
      writer.write(SizeX);
      writer.write("\r\n");
      writer.write(SizeY);
      writer.write("\r\n");

      writer.write(freq);
      writer.write("\r\n");
      writer.write(Nlambda);
      writer.write("\r\n");
      writer.write(maxTime);
      writer.write("\r\n");

      writer.write(TFSF_x0);
      writer.write("\r\n");
      writer.write(TFSF_sizeX);
      writer.write("\r\n");
      writer.write(TFSF_y0);
      writer.write("\r\n");
      writer.write(TFSF_sizeY);
      writer.write("\r\n");
      writer.write(TFSFangle);
      writer.write("\r\n");
      writer.write(srcType);
      writer.write("\r\n");
      writer.write(src);
      writer.write("\r\n");
      writer.write(PML_Size); 
      writer.close();

    } catch (IOException e){
        e.printStackTrace();
    }

    System.exit(0); //Close GUI
  }//end print func
}//end class













