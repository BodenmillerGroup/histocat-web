"""Add member table

Revision ID: ee2d279cae47
Revises: e0e586eec464
Create Date: 2020-09-30 11:19:13.281166

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = 'ee2d279cae47'
down_revision = 'e0e586eec464'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'member',
        sa.Column('id', sa.Integer(), primary_key=True, index=True),
        sa.Column('group_id', sa.Integer(), sa.ForeignKey("group.id", ondelete="CASCADE"), index=True, nullable=False),
        sa.Column('user_id', sa.Integer(), sa.ForeignKey("user.id", ondelete="CASCADE"), index=True, nullable=False),
        sa.Column('role', sa.Integer(), default=0, nullable=False),
        sa.Column('is_active', sa.Boolean(), default=False, nullable=False),
        sa.Column('created_at', sa.DateTime(), default=sa.sql.func.now(), nullable=False),
        sa.Column('updated_at', sa.DateTime(), default=sa.sql.func.now(), onupdate=sa.sql.func.now(), nullable=False),
    )


def downgrade():
    op.drop_table('member')
