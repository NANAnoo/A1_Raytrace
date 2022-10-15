//////////////////////////////////////////////////////////////////////
//
//  University of Leeds
//  COMP 5812M Foundations of Modelling & Rendering
//  User Interface for Coursework
////////////////////////////////////////////////////////////////////////

// system libraries
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

// QT
#include <QApplication>

// local includes
#include "RenderWindow.h"
#include "ThreeDModel.h"
#include "RenderParameters.h"
#include "RenderController.h"

// main routine
int main(int argc, char **argv)
    { // main()
    // initialize QT
    QApplication renderApp(argc, argv);

    // check the args to make sure there's an input file
    if (argc != 3)
        { // bad arg count
        // print an error message
        std::cout << "Usage: " << argv[0] << " geometry texture|material" << std::endl;
        // and leave
        return 0;
        } // bad arg count

    //Our obj file may contain multiple objects.

    std::vector<ThreeDModel> texturedObjects;

    // open the input files for the geometry & texture
    std::ifstream geometryFile(argv[1]);
    std::ifstream textureFile(argv[2]);


    // try reading it
    if (!(geometryFile.good()) || !(textureFile.good())){
        std::cout << "Read failed for object " << argv[1] << " or texture " << argv[2] << std::endl;
        return 0;
        } // object read failed

    std::string s = argv[2];
    //if is actually passing a material. This will trigger the modified obj read code.
    if(s.find(".mtl") != std::string::npos){
        texturedObjects = ThreeDModel::ReadObjectStreamMaterial(geometryFile, textureFile);
    }else{
        texturedObjects = ThreeDModel::ReadObjectStream(geometryFile);
    }

    if(texturedObjects.size() == 0){
        std::cout << "Read failed for object " << argv[1] << " or texture " << argv[2] << std::endl;
        return 0;
    } // object read failed



    // dump the file to out
//      texturedObject.WriteObjectStream(std::cout, std::cout);

    // create some default render parameters
    RenderParameters renderParameters;

    renderParameters.findLights(texturedObjects);

    // use the object & parameters to create a window
    RenderWindow renderWindow(&texturedObjects, &renderParameters, argv[1]);

    // create a controller for the window
    RenderController renderController(&texturedObjects, &renderParameters, &renderWindow);

    //  set the initial size
    renderWindow.resize(1600, 720);

    // show the window
    renderWindow.show();

    // set QT running
    return renderApp.exec();
    } // main()
