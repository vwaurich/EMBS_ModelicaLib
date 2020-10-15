within EMBSlib;
model Beam "A simple beam model"
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
   Modelica.Mechanics.Rotational.Sources.Position position1(useSupport=true)
     annotation (Placement(transformation(extent={{-85,55},{-65,75}})));
   Modelica.Blocks.Sources.RealExpression realExpression1(y=0.5*sin(time))
     annotation (Placement(transformation(extent={{-135,55},{-115,75}})));
 Components.EMBS_Body eMBS_Body1(
    numNodes=8,
    numModes=4,
    SIDfileName=Modelica.Utilities.Files.loadResource(
        "modelica://EMBSlib/Resources/Data/Rues.SID_FEM"))                annotation(Placement(transformation(extent={{-4,0},{
            16,20}})));
 Components.EMBS_Body eMBS_Body2(
    numNodes=8,
    numModes=4,
    SIDfileName=Modelica.Utilities.Files.loadResource(
        "modelica://EMBSlib/Resources/Data/Rues.SID_FEM"))                annotation(Placement(transformation(extent={{48,-8},
            {68,12}})));
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
  connect(eMBS_Body.frame_node[4], eMBS_Body1.frame_ref) annotation (Line(
      points={{-22,9.8},{-14,9.8},{-14,10},{-4,10}},
      color={95,95,95},
      thickness=0.5));
  connect(eMBS_Body1.frame_node[4], eMBS_Body2.frame_ref) annotation (Line(
      points={{16,9.8},{34,9.8},{34,2},{48,2}},
      color={95,95,95},
      thickness=0.5));
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
               100,100}})));
end Beam;
