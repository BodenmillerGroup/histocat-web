from sqlalchemy import inspect
from sqlalchemy.ext.declarative import as_declarative


@as_declarative()
class Base:
    def as_dict(self) -> dict:
        return {c.key: getattr(self, c.key) for c in inspect(self).mapper.column_attrs}
