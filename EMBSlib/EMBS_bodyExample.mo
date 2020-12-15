within EMBSlib;
model EMBS_bodyExample
 extends Modelica.Icons.Example;
 Components.EMBS_Body eMBS_Body(SIDfileName=Modelica.Utilities.Files.loadResource(
  "modelica://EMBSlib/Resources/Data/cartopPragV32.SID_FEM")) annotation(Placement(transformation(extent={{26,-10},
            {46,10}})));
inner Modelica.Mechanics.MultiBody.World world(
 g=9.81,
 n(displayUnit="1")={0,0,-1}) annotation(Placement(transformation(extent={{-54,-10},
            {-34,10}})));
Modelica.Mechanics.MultiBody.Joints.Revolute revFix(
 n(displayUnit="1")={1,
          0,0},
w(
 fixed=true,
 start=0)) annotation(Placement(transformation(extent={{-12,-10},{8,10}})));
equation
 connect(world.frame_b, revFix.frame_a) annotation (Line(
     points={{-34,0},{-12,0}},
     color={95,95,95},
     thickness=0.5));
 connect(eMBS_Body.frame_ref, revFix.frame_b) annotation (Line(
     points={{26,0},{8,0}},
     color={95,95,95},
     thickness=0.5));
annotation(experiment(
 StopTime=10,
 __Dymola_Algorithm="Dassl"));
end EMBS_bodyExample;
