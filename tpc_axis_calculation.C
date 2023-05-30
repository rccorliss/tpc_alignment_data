void tpc_axis_calculation(){
  //this code calculates the transformation needed to convert a magnet map measured with respect to the magnet axes into a field map with respect to the TPC coordinate system.  It hard-codes the survey measurements made of the magnet and TPC by the BNL survey group, and accounts for pitch, yaw, roll, and translation of both the TPC and the magnet with respect to the sPHENIX global coordinate system.

  //These values are used in the AnnularFieldSim code that calculates/simulates distortions in the TPC, and may be useful elsewhere as well.
  
  //we assume our coordinate system is xhat vertical, yhat horizontal and perp to the beam, and zhat along the beam
  //this differs from zurfvey which has zhat vertical, xhat perp to the beam, and y hat along the beam

  //measured in inches
  const float inches=2.54; //inches to cm.
  const float mm=0.10; //mm to cm.
  //magnet position data from Dan Cacace's May 26 email to Ross
  TVector3 magnetPosInGlobal(0.057425*inches,0.079889*inches, 0.177971*inches); //May 26 email data, from magnet survey without TPC
  TVector3 magnetPosInGlobal2(0.638*mm,2.047*mm,4.333*mm); //May 5 email data, from survey with TPC installed
  float magnetPitchInGlobal=-0.0098*TMath::Pi()/180.;
  float magnetRollInGlobal=-0.0004*TMath::Pi()/180.;
  float magnetYawInGlobal=0.0485*TMath::Pi()/180;

  //now derive the magnet axes:
  TVector3 magnetHatInGlobal[3];
  magnetHatInGlobal[0].SetXYZ(1,0,0);
  magnetHatInGlobal[1].SetXYZ(0,1,0);
  magnetHatInGlobal[2].SetXYZ(0,0,1);

  //apply pitch in global:
  magnetHatInGlobal[0].RotateY(magnetPitchInGlobal);
  magnetHatInGlobal[2].RotateY(magnetPitchInGlobal);
  //apply yaw in global:
  magnetHatInGlobal[1].RotateX(magnetYawInGlobal);
  magnetHatInGlobal[2].RotateX(magnetYawInGlobal);
  //apply roll in local:
  magnetHatInGlobal[0].Rotate(magnetRollInGlobal, magnetHatInGlobal[2]);
  magnetHatInGlobal[1].Rotate(magnetRollInGlobal, magnetHatInGlobal[2]);

  //finally, we generate the two rotations to get into and out of magnetHat axes:
  TRotation magnetActiveGlobal; //the transformation that takes global zhat and points it in the magnet zhat direction, ditto xhat.  Equivalently, the transformation that expresses a magnet position and expresses it in terms of the global axes.
  magnetActiveGlobal.SetZAxis(magnetHatInGlobal[2],magnetHatInGlobal[0]);
  TRotation magnetPassiveGlobal; //the transformation that expresses a global position in terms of magnet axis projections.
  magnetPassiveGlobal=magnetActiveGlobal.Inverse();

  //TPC position data from Dan Cacace's May 5 2023 email to Ross
  TVector3 tpcEndInGlobal(-0.673*mm,-3.354*mm,1137.382*mm);
  TVector3 tpcBeginInGlobal(0.001*mm,-0.001*mm,-1123.109*mm);
  TVector3 tpcPosInGlobal=(tpcEndInGlobal+tpcBeginInGlobal)*0.5; //average the position.
  TVector3 tpcAxisInGlobal=tpcEndInGlobal-tpcBeginInGlobal;
TVector3 tpcY0=tpcAxisInGlobal;tpcY0.SetY(0);
  TVector3 tpcX0=tpcAxisInGlobal;tpcX0.SetX(0);
  float tpcPitchInGlobal=tpcAxisInGlobal.Angle(tpcX0); //pitch is the angle between our vector and its shadow in the YZ plane;
  float tpcYawInGlobal=tpcAxisInGlobal.Angle(tpcY0); //yaw is the angle between our vector and its shadow in the XZ plane;
  float tpcRollInGlobal=0.000; //no data on this, so start by assuming it's zero.

  TVector3 tpcHatInGlobal[3];
  tpcHatInGlobal[0].SetXYZ(1,0,0);
  tpcHatInGlobal[1].SetXYZ(0,1,0);
  tpcHatInGlobal[2]=tpcAxisInGlobal.Unit(); //we have the right answer for this already

  //only bother to generate the correct x coord, since we don't use Y.
  tpcHatInGlobal[0].RotateY(tpcPitchInGlobal);
  tpcHatInGlobal[0].Rotate(tpcRollInGlobal,tpcHatInGlobal[2]);

   //finally, we generate the two rotations to get into and out of magnetHat axes:
  TRotation tpcActiveGlobal; //the transformation that takes global zhat and points it in the tpc zhat direction, ditto xhat.
  tpcActiveGlobal.SetZAxis(tpcHatInGlobal[2],tpcHatInGlobal[0]);
  TRotation tpcPassiveGlobal; //the transformation that expresses a global position in terms of tpc axis projections.
  tpcPassiveGlobal=tpcActiveGlobal.Inverse();

  //now we can compose, in order to express the magnet position in TPC local coordinates:
  // tpcCoordinates=netOffset+netRotation*magnetCoordinates
  TVector3 magnetPosInTpc=tpcPassiveGlobal*(magnetPosInGlobal-tpcPosInGlobal); //the position of the magnet center expressed in TPC coordinate system 
  TRotation tpcPassiveMagnet=tpcPassiveGlobal*magnetActiveGlobal;//the transformation that expresses a magnet position in terms of tpc axis projections

  
  printf("Magnet->Global:\n");
  for (int i=0;i<3;i++){
    printf("%f\t%f\t%f\n",magnetActiveGlobal(i,0),magnetActiveGlobal(i,1),magnetActiveGlobal(i,2));
  }
  printf("\n");

 printf("Global->TPC:\n");
  for (int i=0;i<3;i++){
    printf("%f\t%f\t%f\n",tpcPassiveGlobal(i,0),tpcPassiveGlobal(i,1),tpcPassiveGlobal(i,2));
  }
  printf("\n");
  
  
 printf("Magnet->TPC:\n");
  for (int i=0;i<3;i++){
    printf("%f\t%f\t%f\n",tpcPassiveMagnet(i,0),tpcPassiveMagnet(i,1),tpcPassiveMagnet(i,2));
  }
  printf("\n");

  printf("Magnet center in TPC coords:\n");
  magnetPosInTpc.Print();

  //this should work, but doesn't in version 6.26/06
  /*
  Double_t netAngle; TVector3 netAxis;
  tpcPassiveMagnet.GetAngleAxis(netAngle,netAxis);
  printf("Magnet->TPC net angle and axis:\n");
  printf("netAngle=%f\n",netAngle);
  netAxis.Print();
  */
  //report the X Euler Angles, which we can use to rebuild this:

  printf("Magnet->TPC Euler angles:  XPhi=%f, xTheta=%f, XPsi=%f\n",tpcPassiveMagnet.GetXPhi(),tpcPassiveMagnet.GetXTheta(),tpcPassiveMagnet.GetXPsi());
  printf("To generate this transformation:\n  TRotation magnetToTpc; magnetToTpc.SetXEulerAngles(%f,%f,%f)\n",tpcPassiveMagnet.GetXPhi(),tpcPassiveMagnet.GetXTheta(),tpcPassiveMagnet.GetXPsi());

  return;
  
}
