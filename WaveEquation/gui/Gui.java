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
  
  private JTextField SizeXfield, SizeYfield;
  private JTextField dxField, maxTimeField;
  private JTextField FreqField;
  private JTextField CFLfield;
 
  private JRadioButton Harmonic, Gaussian;
  private ButtonGroup Source;
  
  private GridBagConstraints c; 
  static public JPanel panel;  
  private JButton doneButton; 
  
  private String SizeX, SizeY;
  private String dx, maxTime;
  private String freq;
  private String CFL, src;
  
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
    placeField("Size X (m) :    ", SizeXfield);
    
    SizeYfield = new JTextField("200     ");
    placeField("Size Y (m) :    ", SizeYfield);
    
    dxField = new JTextField("1.0     ");
    placeField("dx (mm) :    ", dxField);
    
    CFLfield = new JTextField("0.99      ");
    placeField("CFL (stability) :   ", CFLfield);
    
    maxTimeField = new JTextField("100         ");  
    placeField("Max Time (s) : ", maxTimeField);

    FreqField = new JTextField("10        ");  
    placeField("Frequency (GHz) : ", FreqField); 
    
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
      
      dx  = dxField.getText();
      CFL = CFLfield.getText(); 
      maxTime = maxTimeField.getText();

      freq = FreqField.getText();
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
 //   Harmonic = new JRadioButton("Harmonic");
  //  Gaussian = new JRadioButton("Gaussian");

   // Source = new ButtonGroup();
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

  public void printToFile(){
    
    try{
      FileWriter writer = new FileWriter("grid.txt");
      writer.write(SizeX);
      writer.write("\r\n");
      writer.write(SizeY);
      writer.write("\r\n");

      writer.write(dx);
      writer.write("\r\n");
      writer.write(CFL);
      writer.write("\r\n");

      writer.write(maxTime);
      writer.write("\r\n"); 
      writer.write(freq);
      writer.write("\r\n");
      writer.write(src);
      writer.close();

    } catch (IOException e){
        e.printStackTrace();
    }

    System.exit(0); //Close GUI
  }//end print func
}//end class













