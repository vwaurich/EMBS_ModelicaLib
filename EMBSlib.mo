within ;
package EMBSlib

  class SID_Data
    extends ExternalObject;

    function constructor
      input String fileName;
      output SID_Data sid;
      external "C" sid = SID_Constructor(fileName)
      annotation(Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});

    end constructor;

    function destructor
      input SID_Data sid;
      external "C" SID_Destructor(sid)
      annotation(Include="#include \"SID_Data.h\"",
                   Library={"readSIDlib"});
    end destructor;

  end SID_Data;

  package ExternalFunctions
    function getNumberOfNodes
      input EMBSlib.SID_Data sid;
      output Integer numNodes;
      external "C" numNodes = getNumberOfNodes(sid)
      annotation(Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
    end getNumberOfNodes;

    function getNumberOfModes
      input EMBSlib.SID_Data sid;
      output Integer numNodes;
      external "C" numNodes = getNumberOfModes(sid)
      annotation(Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
    end getNumberOfModes;

    function getOriginOrder
      input EMBSlib.SID_Data sid;
      input Integer nodeIdx;
      output Integer order;
    external"C" order = getOriginOrder(sid, nodeIdx) annotation (
          Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
    end getOriginOrder;

      function getOriginNRows
      input EMBSlib.SID_Data sid;
      input Integer nodeIdx;
      output Integer numRows;
      external "C" numRows = getOriginNRows(sid, nodeIdx)
      annotation(Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
      end getOriginNRows;

    function getOriginNCols
      input EMBSlib.SID_Data sid;
      input Integer nodeIdx;
      output Integer numCols;
    external"C" numCols = getOriginNCols(sid, nodeIdx) annotation (
          Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
    end getOriginNCols;

    function getOriginNq
      input EMBSlib.SID_Data sid;
      input Integer nodeIdx;
      output Integer order;
    external"C" order = getOriginNq(sid, nodeIdx) annotation (
          Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
    end getOriginNq;

      function getOriginNqn
      input EMBSlib.SID_Data sid;
      input Integer nodeIdx;
      output Integer order;
    external"C" order = getOriginNqn(sid, nodeIdx) annotation (
          Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
      end getOriginNqn;

     function getOriginStructure
      input EMBSlib.SID_Data sid;
      input Integer nodeIdx;
      output Integer order;
    external"C" order = getOriginStructure(sid, nodeIdx) annotation (
          Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
     end getOriginStructure;

      function getM0ForNode
      input EMBSlib.SID_Data sid;
      input Integer nodeIdx;
      input Integer numRows;
      input Integer numCols;
      output Real[numRows, numCols] M0;
      external "C" getM0ArrforNode(sid, nodeIdx, M0)
      annotation(Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
      end getM0ForNode;

    function getOriginM0
      input EMBSlib.SID_Data sid;
      input Integer nodeIdx;
      input Integer r;
      input Integer c;
      input Integer dimR;
      input Integer dimC;
      output Real value;
      external "C" value = getOriginM0(sid,nodeIdx,r,c, dimR, dimC)
      annotation(Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
    end getOriginM0;

      function getOriginM1
      input EMBSlib.SID_Data sid;
      input Integer nodeIdx;
      input Integer r;
      input Integer q;
      input Integer c;
      input Integer dimR;
      input Integer dimQ;
      input Integer dimC;
      output Real value;
      external "C" value = getOriginM1(sid,nodeIdx,r,q,c, dimR,dimQ, dimC)
      annotation(Include="#include \"SID_Data.h\"",
                  Library={"readSIDlib"});
      end getOriginM1;
  end ExternalFunctions;
  annotation (uses(Modelica(version="3.2.1")));
  package Types
    record Taylor
      Integer order;
      Integer nrow;
      Integer ncol;
      Integer nq;
      Integer nqn;
      Integer structure;
      Real[3,1] M0; //nrows,ncol
      Real[3,3,1] M1; //rows,q,cols
    end Taylor;

    record Matrix31
        Real M[3];
    end Matrix31;
  end Types;

  package Functions
    function getOrigin
        input EMBSlib.SID_Data sid;
        input Integer nodeIdx;
        output EMBSlib.Types.Taylor nodeOrig;
    algorithm
      nodeOrig.order := EMBSlib.ExternalFunctions.getOriginOrder(sid, nodeIdx);
      nodeOrig.nrow := EMBSlib.ExternalFunctions.getOriginNRows(sid, nodeIdx);
      nodeOrig.ncol := EMBSlib.ExternalFunctions.getOriginNCols(sid, nodeIdx);
      nodeOrig.nq := EMBSlib.ExternalFunctions.getOriginNq(sid, nodeIdx);
      nodeOrig.nqn := EMBSlib.ExternalFunctions.getOriginNqn(sid, nodeIdx);
      nodeOrig.structure := EMBSlib.ExternalFunctions.getOriginStructure(sid,
        nodeIdx);
      nodeOrig.M0[1, 1] := EMBSlib.ExternalFunctions.getOriginM0(
            sid,
            nodeIdx,
            1,
            1,
            nodeOrig.nrow,
            nodeOrig.ncol);
      nodeOrig.M0[2, 1] := EMBSlib.ExternalFunctions.getOriginM0(
            sid,
            nodeIdx,
            2,
            1,
            nodeOrig.nrow,
            nodeOrig.ncol);
      nodeOrig.M0[3, 1] := EMBSlib.ExternalFunctions.getOriginM0(
            sid,
            nodeIdx,
            3,
            1,
            nodeOrig.nrow,
            nodeOrig.ncol);
      nodeOrig.M1[1, 1, 1] := EMBSlib.ExternalFunctions.getOriginM1(
            sid,
            nodeIdx,
            1,
            1,
            1,
            nodeOrig.nrow,
            nodeOrig.nq,
            nodeOrig.ncol);
      nodeOrig.M1[2, 1, 1] := EMBSlib.ExternalFunctions.getOriginM1(
            sid,
            nodeIdx,
            2,
            1,
            1,
            nodeOrig.nrow,
            nodeOrig.nq,
            nodeOrig.ncol);
      nodeOrig.M1[3, 1, 1] := EMBSlib.ExternalFunctions.getOriginM1(
            sid,
            nodeIdx,
            3,
            1,
            1,
            nodeOrig.nrow,
            nodeOrig.nq,
            nodeOrig.ncol);
      nodeOrig.M1[1, 2, 1] := EMBSlib.ExternalFunctions.getOriginM1(
            sid,
            nodeIdx,
            1,
            2,
            1,
            nodeOrig.nrow,
            nodeOrig.nq,
            nodeOrig.ncol);
      nodeOrig.M1[2, 2, 1] := EMBSlib.ExternalFunctions.getOriginM1(
            sid,
            nodeIdx,
            2,
            2,
            1,
            nodeOrig.nrow,
            nodeOrig.nq,
            nodeOrig.ncol);
      nodeOrig.M1[3, 2, 1] := EMBSlib.ExternalFunctions.getOriginM1(
            sid,
            nodeIdx,
            3,
            2,
            1,
            nodeOrig.nrow,
            nodeOrig.nq,
            nodeOrig.ncol);
      nodeOrig.M1[1, 3, 1] := EMBSlib.ExternalFunctions.getOriginM1(
            sid,
            nodeIdx,
            1,
            3,
            1,
            nodeOrig.nrow,
            nodeOrig.nq,
            nodeOrig.ncol);
      nodeOrig.M1[2, 3, 1] := EMBSlib.ExternalFunctions.getOriginM1(
            sid,
            nodeIdx,
            2,
            3,
            1,
            nodeOrig.nrow,
            nodeOrig.nq,
            nodeOrig.ncol);
      nodeOrig.M1[3, 3, 1] := EMBSlib.ExternalFunctions.getOriginM1(
            sid,
            nodeIdx,
            3,
            3,
            1,
            nodeOrig.nrow,
            nodeOrig.nq,
            nodeOrig.ncol);
    end getOrigin;
  end Functions;

  package Components
    model Nodes
      parameter Integer numNodes = 9;
      parameter String SIDfileName = "E:/Projekte/VIBROSIM_2/EMBS_ModelicaLib/Resources/Data/cartopPragV32.SID_FEM";
      EMBSlib.Types.Taylor nodeOrig[numNodes];
      EMBSlib.SID_Data sid = EMBSlib.SID_Data(SIDfileName);

    algorithm
      for i in 1:numNodes loop
          nodeOrig[i].order := EMBSlib.ExternalFunctions.getOriginOrder(sid, i);
      nodeOrig[i].nrow := EMBSlib.ExternalFunctions.getOriginNRows(sid, i);
      nodeOrig[i].ncol := EMBSlib.ExternalFunctions.getOriginNCols(sid, i);
      nodeOrig[i].nq := EMBSlib.ExternalFunctions.getOriginNq(sid, i);
      nodeOrig[i].nqn := EMBSlib.ExternalFunctions.getOriginNqn(sid, i);
      nodeOrig[i].structure := EMBSlib.ExternalFunctions.getOriginStructure(sid,
        i);
      nodeOrig[i].M0[1, 1] := EMBSlib.ExternalFunctions.getOriginM0(
        sid,
        i,
        1,
        1,
        nodeOrig[i].nrow,
        nodeOrig[i].ncol);
      nodeOrig[i].M0[2, 1] := EMBSlib.ExternalFunctions.getOriginM0(
        sid,
        i,
        2,
        1,
        nodeOrig[i].nrow,
        nodeOrig[i].ncol);
      nodeOrig[i].M0[3, 1] := EMBSlib.ExternalFunctions.getOriginM0(
        sid,
        i,
        3,
        1,
        nodeOrig[i].nrow,
        nodeOrig[i].ncol);
      nodeOrig[i].M1[1, 1, 1] := EMBSlib.ExternalFunctions.getOriginM1(
        sid,
        i,
        1,
        1,
        1,
        nodeOrig[i].nrow,
        nodeOrig[i].nq,
        nodeOrig[i].ncol);
      nodeOrig[i].M1[2, 1, 1] := EMBSlib.ExternalFunctions.getOriginM1(
        sid,
        i,
        2,
        1,
        1,
        nodeOrig[i].nrow,
        nodeOrig[i].nq,
        nodeOrig[i].ncol);
      nodeOrig[i].M1[3, 1, 1] := EMBSlib.ExternalFunctions.getOriginM1(
        sid,
        i,
        3,
        1,
        1,
        nodeOrig[i].nrow,
        nodeOrig[i].nq,
        nodeOrig[i].ncol);
      nodeOrig[i].M1[1, 2, 1] := EMBSlib.ExternalFunctions.getOriginM1(
        sid,
        i,
        1,
        2,
        1,
        nodeOrig[i].nrow,
        nodeOrig[i].nq,
        nodeOrig[i].ncol);
      nodeOrig[i].M1[2, 2, 1] := EMBSlib.ExternalFunctions.getOriginM1(
        sid,
        i,
        2,
        2,
        1,
        nodeOrig[i].nrow,
        nodeOrig[i].nq,
        nodeOrig[i].ncol);
      nodeOrig[i].M1[3, 2, 1] := EMBSlib.ExternalFunctions.getOriginM1(
        sid,
        i,
        3,
        2,
        1,
        nodeOrig[i].nrow,
        nodeOrig[i].nq,
        nodeOrig[i].ncol);
      nodeOrig[i].M1[1, 3, 1] := EMBSlib.ExternalFunctions.getOriginM1(
        sid,
        i,
        1,
        3,
        1,
        nodeOrig[i].nrow,
        nodeOrig[i].nq,
        nodeOrig[i].ncol);
      nodeOrig[i].M1[2, 3, 1] := EMBSlib.ExternalFunctions.getOriginM1(
        sid,
        i,
        2,
        3,
        1,
        nodeOrig[i].nrow,
        nodeOrig[i].nq,
        nodeOrig[i].ncol);
      nodeOrig[i].M1[3, 3, 1] := EMBSlib.ExternalFunctions.getOriginM1(
        sid,
        i,
        3,
        3,
        1,
        nodeOrig[i].nrow,
        nodeOrig[i].nq,
        nodeOrig[i].ncol);
      end for;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
                100,100}}), graphics));

    end Nodes;

    model EMBS_Body
      parameter Integer numNodes = 9;

      Nodes nodes(numNodes=numNodes)
                  annotation (Placement(transformation(extent={{0,40},{20,60}})));

    Real              rShapes[numNodes,3];
      inner Modelica.Mechanics.MultiBody.World world
        annotation (Placement(transformation(extent={{-100,40},{-80,60}})));
      Modelica.Mechanics.MultiBody.Visualizers.Advanced.Shape shape[numNodes](r= rShapes,
        shapeType="sphere",
        length=0.1,
        width=0.1,
        height=0.1)
        annotation (Placement(transformation(extent={{-60,-20},{-40,0}})));
    algorithm
        for i in 1:numNodes loop
          rShapes[i,1] := nodes.nodeOrig[i].M0[1,1];
          rShapes[i,2] := nodes.nodeOrig[i].M0[2,1];
          rShapes[i,3] := nodes.nodeOrig[i].M0[3,1];
        end for;
    end EMBS_Body;
  end Components;
end EMBSlib;
