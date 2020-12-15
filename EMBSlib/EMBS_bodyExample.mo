within EMBSlib;
model EMBS_bodyExample
 extends Modelica.Icons.Example;
 Components.EMBS_Body eMBS_Body(
    numNodes=6,
    numModes=11,
    SIDfileName=Modelica.Utilities.Files.loadResource(
        "modelica://EMBSlib/Resources/Data/beam.SID_FEM"))    annotation(Placement(transformation(extent={{0,0},{
            20,20}})));
inner Modelica.Mechanics.MultiBody.World world(
    label2="y",
 g=9.81,
    n(displayUnit="1") = {0,-1,0})
                              annotation(Placement(transformation(extent={{-40,0},
            {-20,20}})));
  Modelica.Mechanics.MultiBody.Forces.WorldForce force annotation (Placement(
        transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={50,30})));
  Modelica.Blocks.Sources.RealExpression zero
    annotation (Placement(transformation(extent={{20,60},{40,80}})));
  Modelica.Blocks.Sources.RealExpression zSignal(y=-time*10000)
    annotation (Placement(transformation(extent={{20,40},{40,60}})));
equation
  connect(eMBS_Body.frame_ref, world.frame_b) annotation (Line(
      points={{0,10},{-20,10}},
      color={95,95,95},
      thickness=0.5));
  connect(force.frame_b, eMBS_Body.frame_node[3]) annotation (Line(
      points={{50,20},{50,10},{36,10},{36,9.73333},{20,9.73333}},
      color={95,95,95},
      thickness=0.5));
  connect(zero.y, force.force[1]) annotation (Line(points={{41,70},{60,70},{60,
          42},{48.6667,42}}, color={0,0,127}));
  connect(zSignal.y, force.force[2])
    annotation (Line(points={{41,50},{50,50},{50,42}}, color={0,0,127}));
  connect(zero.y, force.force[3]) annotation (Line(points={{41,70},{60,70},{60,
          42},{51.3333,42}}, color={0,0,127}));
annotation(experiment(StopTime=10, __Dymola_Algorithm="Dassl"));
end EMBS_bodyExample;
