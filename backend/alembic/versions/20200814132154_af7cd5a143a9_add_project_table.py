"""Add project table

Revision ID: af7cd5a143a9
Revises: e1c9291e68ae
Create Date: 2020-08-14 13:21:54.388521

"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import ARRAY


# revision identifiers, used by Alembic.
revision = 'af7cd5a143a9'
down_revision = 'e1c9291e68ae'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'project',
        sa.Column('id', sa.Integer(), primary_key=True, index=True),
        sa.Column('group_id', sa.Integer(), sa.ForeignKey("group.id", ondelete="CASCADE"), index=True, nullable=False),
        sa.Column('member_id', sa.Integer(), sa.ForeignKey("member.id", ondelete="CASCADE"), index=True, nullable=False),
        sa.Column('name', sa.String(), nullable=False, index=True, unique=True),
        sa.Column('description', sa.String()),
        sa.Column('tags', ARRAY(sa.String(64)), index=True),
        sa.Column('location', sa.String(4096)),
        sa.Column('created_at', sa.DateTime(), default=sa.sql.func.now(), nullable=False),
    )


def downgrade():
    op.drop_table('project')
