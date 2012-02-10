
ADD_EXECUTABLE( test1d Test1d.cxx )
TARGET_LINK_LIBRARIES ( test1d
vul 
${ITK_LIBRARIES})

ADD_EXECUTABLE( rpcfm RegisterParameterizedCFM.cxx ParameterInitialization.cxx ParameterizedLevelset.cxx)
TARGET_LINK_LIBRARIES ( rpcfm
vul 
json
${ITK_LIBRARIES})
