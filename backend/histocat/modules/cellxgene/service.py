import copy
import logging

logger = logging.getLogger(__name__)


def get_schema(data_adaptor, annotations):
    """helper function to gather the schema from the data source and annotations"""
    schema = data_adaptor.get_schema()
    schema = copy.deepcopy(schema)

    # add label obs annotations as needed
    if annotations is not None:
        label_schema = annotations.get_schema(data_adaptor)
        schema["annotations"]["obs"]["columns"].extend(label_schema)

    return schema
