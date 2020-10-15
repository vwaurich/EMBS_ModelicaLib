within EMBSlib;
model EMBS_bodyExample
 extends Modelica.Icons.Example;
 Components.EMBS_Body eMBS_Body(SIDfileName=Modelica.Utilities.Files.loadResource(
  "modelica://EMBSlib/Resources/Data/cartopPragV32.SID_FEM")) annotation(Placement(transformation(extent={{-20,0},{0,20}})));
inner Modelica.Mechanics.MultiBody.World world(
 g=9.81,
 n(displayUnit="1")={0,0,-1}) annotation(Placement(transformation(extent={{-100,0},{-80,20}})));
Modelica.Mechanics.MultiBody.Joints.Revolute revFix(
 n(displayUnit="1")={1,
          0,0},
w(
 fixed=true,
 start=0)) annotation(Placement(transformation(extent={{-58,0},{-38,20}})));
equation
 connect(world.frame_b, revFix.frame_a) annotation (Line(
     points={{-80,10},{-58,10}},
     color={95,95,95},
     thickness=0.5));
 connect(eMBS_Body.frame_ref, revFix.frame_b) annotation (Line(
     points={{-20,10},{-38,10}},
     color={95,95,95},
     thickness=0.5));
annotation(experiment(
 StopTime=10,
 __Dymola_Algorithm="Dassl"));
end EMBS_bodyExample;
