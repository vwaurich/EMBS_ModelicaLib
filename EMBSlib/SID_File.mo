within EMBSlib;
class SID_File
  extends ExternalObject;

  function constructor
    input String fileName;
    output SID_File sid;
    external "C" sid = SIDFileConstructor_C(fileName)
    annotation(Include="#include \"ReadSID_C.h\"");

  end constructor;

  function destructor
    input SID_File sid;
    external "C" SIDFileDestructor_C(sid)
    annotation(Include="#include \"ReadSID_C.h\"");
  end destructor;

end SID_File;
