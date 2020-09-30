"""Add preset table

Revision ID: 59a8939464a4
Revises: b3c290b4bb05
Create Date: 2020-09-30 11:22:41.822038

"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import JSONB


# revision identifiers, used by Alembic.
revision = '59a8939464a4'
down_revision = 'b3c290b4bb05'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'preset',
        sa.Column('id', sa.Integer(), primary_key=True, index=True),
        sa.Column('project_id', sa.Integer(), sa.ForeignKey("project.id", ondelete="CASCADE"), index=True, nullable=False),
        sa.Column('member_id', sa.Integer(), sa.ForeignKey("member.id", ondelete="CASCADE"), index=True, nullable=False),
        sa.Column('name', sa.String()),
        sa.Column('description', sa.String()),
        sa.Column('data', JSONB()),
        sa.Column('created_at', sa.DateTime(), default=sa.sql.func.now(), nullable=False),
    )


def downgrade():
    op.drop_table('preset')
