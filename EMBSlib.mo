within ;
package EMBSlib
  model testSID
    Integer i1;
    EMBSlib.SID_Data sid = EMBSlib.SID_Data("aha");
  algorithm
     i1 :=EMBSlib.Functions.getNumberOfNodes(sid);
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end testSID;

  class SID_Data
    extends ExternalObject;

    function constructor
      input String fileName;
      output SID_Data sid;
    external "C" sid = SID_Constructor(fileName) annotation(Include="#include \"SID_Data.h\"");
    end constructor;

    function destructor
      input SID_Data sid;
    external "C" SID_Destructor(sid) annotation(Include="#include \"SID_Data.h\"");
    end destructor;

  end SID_Data;

  package Functions
    function getNumberOfNodes
      input EMBSlib.SID_Data sid;
      output Integer numNodes;
      external "C" numNodes = getNumberOfNodes(sid) annotation(Include="#include \"SID_Data.h\"");
    end getNumberOfNodes;
  end Functions;
end EMBSlib;
