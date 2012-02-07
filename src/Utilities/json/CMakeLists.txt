
SET(json_sources
autolink.h
config.h
features.h
forwards.h
json.h
json_batchallocator.h
json_internalarray.inl
json_internalmap.inl
json_reader.cpp
json_value.cpp
json_valueiterator.inl
json_writer.cpp
json_reader.h
json_value.h
json_writer.h
)

ADD_LIBRARY( json ${json_sources} )

