package org.bigbio.pgatk.pepgenome.gui;

import javafx.application.Application;
import javafx.fxml.FXMLLoader;
import javafx.scene.Parent;
import javafx.scene.Scene;
import javafx.scene.layout.VBox;
import javafx.stage.Stage;

public class PepGenomeGUI extends Application {


    @Override
    public void start(Stage stage) throws Exception {

        VBox root = (VBox) FXMLLoader.load(getClass().getClassLoader().getResource("view/PepGenome.fxml"));
        Scene scene = new Scene(root, 640, 532);

        stage.setTitle("PepGenome Tool");
        stage.setScene(scene);
        stage.show();

    }

    public static void main(String[] args){
        Application.launch(args);
    }


}
