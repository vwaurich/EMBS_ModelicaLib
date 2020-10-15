within EMBSlib;
model Balken "Balken"
 extends Modelica.Icons.Example;
 Components.EMBS_Body eMBS_Body(
  numNodes=8,
  numModes=4,
  SIDfileName=Modelica.Utilities.Files.loadResource(
                       "modelica://EMBSlib/Resources/Data/Rues.SID_FEM")) annotation(Placement(transformation(extent={{-42,0},{-22,20}})));
inner Modelica.Mechanics.MultiBody.World world(
 label2="z",
 g=9.81,
 n(displayUnit="1")={0,0,-1}) annotation(Placement(transformation(extent={{-120,0},{-100,20}})));
Modelica.Mechanics.MultiBody.Joints.Revolute revolute1(
 useAxisFlange=true,
 n={0,1,0},
 w(fixed=false)) annotation(Placement(transformation(extent={{-78,0},{-58,20}})));
Components.EMBS_Body eMBS_Body1(
 numNodes=8,
 numModes=4,
 SIDfileName=Modelica.Utilities.Files.loadResource(
                        "modelica://EMBSlib/Resources/Data/Rues.SID_FEM")) annotation(Placement(transformation(extent={{35,-5},{55,15}})));
   Modelica.Mechanics.Rotational.Sources.Position position1(useSupport=true)
     annotation (Placement(transformation(extent={{-85,55},{-65,75}})));
   Modelica.Blocks.Sources.RealExpression realExpression1(y=time*3)
     annotation (Placement(transformation(extent={{-135,55},{-115,75}})));
   Modelica.Mechanics.MultiBody.Joints.Revolute revolute2(
useAxisFlange=true,
n={0,1,0}) annotation(Placement(transformation(extent={{0,0},{20,20}})));
   Modelica.Mechanics.Rotational.Sources.Position position2(useSupport=true)
     annotation (Placement(transformation(extent={{-5,55},{15,75}})));
equation
connect(revolute1.frame_a,world.frame_b) annotation(Line(
 points={{-78,10},{-78,10},{-95,10},{-100,10}},
 color={95,95,95},
 thickness=0.0625));
connect(revolute1.frame_b,eMBS_Body.frame_ref) annotation(Line(
 points={{-58,10},{-58,10},{-47,10},{-42,10}},
 color={95,95,95},
 thickness=0.0625));
connect(position1.support,revolute1.support) annotation(Line(
 points={{-75,55},{-75,50},{-75,25},{-74,25},{-74,20}},
 thickness=0.0625));
connect(position1.flange,revolute1.axis) annotation(Line(
 points={{-65,65},{-60,65},{-60,45},{-70,45},{-70,20},{-68,20}},
 thickness=0.0625));
connect(realExpression1.y,position1.phi_ref) annotation(Line(
 points={{-114,65},{-109,65},{-92,65},{-87,65}},
 color={0,0,127},
 thickness=0.0625));
connect(position2.flange,revolute2.axis) annotation(Line(
 points={{15,65},{20,65},{20,45},{10,45},{10,25},{10,
 20}},
 thickness=0.0625));
connect(position2.support,revolute2.support) annotation(Line(
 points={{5,55},{5,50},{5,25},{4,25},{4,20}},
 thickness=0.0625));
connect(eMBS_Body.frame_node[4],revolute2.frame_a) annotation(Line(
 points={{-22,9.8},{-17,10},{-5,10},{0,10}},
 color={95,95,95},
 thickness=0.0625));
connect(revolute2.frame_b,eMBS_Body1.frame_ref) annotation(Line(
 points={{20,10},{25,10},{30,10},{30,5},{35,5}},
 color={95,95,95},
 thickness=0.0625));
connect(position2.phi_ref,realExpression1.y) annotation(Line(
 points={{-7,65},{-12,65},{-114,65}},
 color={0,0,127},
 thickness=0.0625));
   annotation (
experiment(
 StopTime=1,
 StartTime=0,
 Interval=0.002,
 Algorithm="Dassl",
 Solver(
  animRecord=true,
  typename="BDFCompiled"),
 SolverOptions(
  solver="BDFCompiled",
  typename="BDFCompiledOptionData")),
Diagram(
 coordinateSystem(
  preserveAspectRatio=false,
  extent={{-100,-100},{
               100,100}}),
graphics));
end Balken;
