"""Add member table

Revision ID: e1c9291e68ae
Revises: 4bec27c370f0
Create Date: 2020-08-14 13:21:46.721573

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = 'e1c9291e68ae'
down_revision = '4bec27c370f0'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'member',
        sa.Column('id', sa.Integer(), primary_key=True, index=True),
        sa.Column('group_id', sa.Integer(), sa.ForeignKey("group.id", ondelete="CASCADE"), index=True),
        sa.Column('user_id', sa.Integer(), sa.ForeignKey("user.id", ondelete="CASCADE"), index=True),
        sa.Column('role', sa.Integer(), default=0, nullable=False),
        sa.Column('is_active', sa.Boolean(), default=False, nullable=False),
        sa.Column('created_at', sa.DateTime(), default=sa.sql.func.now(), nullable=False),
        sa.Column('updated_at', sa.DateTime(), default=sa.sql.func.now(), onupdate=sa.sql.func.now(), nullable=False),
    )


def downgrade():
    op.drop_table('member')
