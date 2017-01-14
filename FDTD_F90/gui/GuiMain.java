import javax.swing.*;

public class GuiMain {
  public static void main(String[] args) {
    JFrame frame = new JFrame("FDTD");
    Gui foo = new Gui();   
    frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    frame.getContentPane().add(foo.panel);
    frame.setSize(300,400); 
    frame.setVisible(true);
  }
}//end class
