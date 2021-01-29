"""Add project table

Revision ID: b33b3a0d8331
Revises: ee2d279cae47
Create Date: 2020-09-30 11:19:50.982825

"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import ARRAY


# revision identifiers, used by Alembic.
revision = 'b33b3a0d8331'
down_revision = 'ee2d279cae47'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'project',
        sa.Column('id', sa.Integer(), primary_key=True, index=True),
        sa.Column('group_id', sa.Integer(), sa.ForeignKey("group.id", ondelete="CASCADE"), index=True, nullable=False),
        sa.Column('member_id', sa.Integer(), sa.ForeignKey("member.id", ondelete="CASCADE"), index=True, nullable=False),
        sa.Column('name', sa.String(), nullable=False, index=True),
        sa.Column('description', sa.String()),
        sa.Column('tags', ARRAY(sa.String(64)), index=True),
        sa.Column('location', sa.String(4096)),
        sa.Column('created_at', sa.DateTime(), default=sa.sql.func.now(), nullable=False),
        sa.UniqueConstraint('group_id', 'name', name='uix_project_group_id_and_name')
    )


def downgrade():
    op.drop_table('project')
