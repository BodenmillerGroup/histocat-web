import sqlalchemy as sa

from histocat.db.base import Base


class MemberModel(Base):
    __tablename__ = "member"

    id = sa.Column(sa.Integer(), primary_key=True, index=True)
    group_id = sa.Column(sa.Integer(), sa.ForeignKey("group.id", ondelete="CASCADE"), index=True)
    user_id = sa.Column(sa.Integer(), sa.ForeignKey("user.id", ondelete="CASCADE"), index=True)
    role = sa.Column(sa.Integer(), default=0, nullable=False)
    is_active = sa.Column(sa.Boolean(), default=False, nullable=False)
    created_at = sa.Column(sa.DateTime(), default=sa.sql.func.now(), nullable=False)
    updated_at = sa.Column(sa.DateTime(), default=sa.sql.func.now(), onupdate=sa.sql.func.now(), nullable=False)

    group = sa.orm.relationship("GroupModel", back_populates="members")
    user = sa.orm.relationship("UserModel", back_populates="members")
    projects = sa.orm.relationship("ProjectModel", back_populates="member", cascade="all, delete, delete-orphan")
    presets = sa.orm.relationship("PresetModel", back_populates="member", cascade="all, delete, delete-orphan")

    def __str__(self):
        return f"<{self.__class__.__name__}(id={self.id}, group_id={self.group_id}, , user_id={self.user_id})>"
