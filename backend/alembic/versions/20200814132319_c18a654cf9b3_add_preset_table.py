"""Add preset table

Revision ID: c18a654cf9b3
Revises: a01d5224936c
Create Date: 2020-08-14 13:23:19.279740

"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import JSONB


# revision identifiers, used by Alembic.
revision = 'c18a654cf9b3'
down_revision = 'a01d5224936c'
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
