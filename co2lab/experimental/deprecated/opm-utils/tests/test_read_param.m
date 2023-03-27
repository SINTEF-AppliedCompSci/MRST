%s=readParam(['../data/case.param']);
s=readParam(['../data/spu_2p.param'])% test output 
s=readParam(['../data/case.param'])% test input normal
a=readXmlParam(['../data/case.xml'])% test input xml
s=paramToStruct(['../data/spu_2p.param'])
paramStructToParamFile(s,'spu_write.para')
s=paramToStruct(['../data/case.param'])
a=paramXmlToStruct(['../data/case.xml'])
paramStructToParamFile(a,'case_xml_write.para')
