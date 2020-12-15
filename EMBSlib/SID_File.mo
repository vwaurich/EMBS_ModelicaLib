within EMBSlib;
class SID_File
 extends ExternalObject;
 function constructor
  input String fileName;
  output SID_File sid;

  external "C" sid=SIDFileConstructor_C(fileName) annotation(Include="#include \"ReadSID_C.h\"");
 end constructor;

 function destructor
  input SID_File sid;

  external "C" SIDFileDestructor_C(sid) annotation(Include="#include \"ReadSID_C.h\"");
 end destructor;
  annotation (Icon(graphics={
        Line(points={{-60,-80},{60,-80},{60,40},{20,80},{-60,80},{-60,-80}},
            color={0,0,0}),
        Line(points={{20,80},{20,40},{60,40}}, color={0,0,0}),
        Text(
          extent={{-48,-4},{48,-42}},
          lineColor={0,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.None,
          textString=".SID_FEM")}));
end SID_File;
