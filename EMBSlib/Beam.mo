within EMBSlib;
model Beam "A simple beam model"
 extends Modelica.Icons.Example;
 Components.EMBS_Body eMBS_Body(
    numNodes=3,
    numModes=1,
    SIDfileName=Modelica.Utilities.Files.loadResource(
        "modelica://EMBSlib/Resources/Data/Balken1m.SID_FEM"))            annotation(Placement(transformation(extent={{0,0},{
            20,20}})));
inner Modelica.Mechanics.MultiBody.World world(
 label2="z",
 g=9.81,
 n(displayUnit="1")={0,0,-1}) annotation(Placement(transformation(extent={{-60,0},
            {-40,20}})));
  Modelica.Mechanics.MultiBody.Parts.BodyShape bodyShape(
    r={0.1,0,0},
    r_CM={0.1,0,0},
    m=5,
    length=0.05,
    width=0.05,
    height=0.05) annotation (Placement(transformation(extent={{40,0},{60,20}})));
equation
  connect(bodyShape.frame_a, eMBS_Body.frame_node[1]) annotation (Line(
      points={{40,10},{40,8.93333},{18.3333,8.93333}},
      color={95,95,95},
      thickness=0.5));
  connect(eMBS_Body.frame_ref, world.frame_b) annotation (Line(
      points={{1.66667,10},{-40,10}},
      color={95,95,95},
      thickness=0.5));
   annotation (
experiment(
      StopTime=10,
      Interval=0.002,
      __Dymola_Algorithm="Dassl"),
Diagram(
 coordinateSystem(
  preserveAspectRatio=false,
  extent={{-100,-100},{
               100,100}})));
end Beam;
