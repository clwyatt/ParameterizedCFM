# Testing/CMakeLists.txt

#######################################################
# 2D Test using Digital Phantom with topological change
#######################################################

# build the testfile generator
ADD_EXECUTABLE( make-test-images-2d make-test-images-2d.cxx )
TARGET_LINK_LIBRARIES( make-test-images-2d
${ITK_LIBRARIES} vul)

# setup generator command
ADD_CUSTOM_COMMAND(
OUTPUT circle.nrrd circle-hole.nrrd
COMMAND make-test-images-2d ${CMAKE_CURRENT_BINARY_DIR}/circle.nrrd 
                            ${CMAKE_CURRENT_BINARY_DIR}/circle-hole.nrrd 
DEPENDS make-test-images-2d
)

# trigger the generator
ADD_CUSTOM_TARGET(make-test-images-2d-driver ALL DEPENDS
${CMAKE_CURRENT_BINARY_DIR}/circle.nrrd
${CMAKE_CURRENT_BINARY_DIR}/circle-hole.nrrd
)

# do the test
ADD_TEST(test_images_2d_simple
         ${EXECUTABLE_OUTPUT_PATH}/rpcfm
         ${CMAKE_CURRENT_BINARY_DIR}/circle-hole.nrrd
         ${CMAKE_CURRENT_BINARY_DIR}/circle.nrrd
	 -c ${CMAKE_CURRENT_SOURCE_DIR}/test_images_2d_simple.json
         -o ${CMAKE_CURRENT_BINARY_DIR}/output.nrrd
)
