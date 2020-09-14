"""Add dataset table

Revision ID: a01d5224936c
Revises: 601250049059
Create Date: 2020-08-14 13:22:42.519548

"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy import text
from sqlalchemy.dialects.postgresql import UUID, JSONB


# revision identifiers, used by Alembic.
revision = 'a01d5224936c'
down_revision = '601250049059'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'result',
        sa.Column('id', sa.Integer(), primary_key=True, index=True),
        sa.Column('dataset_id', sa.Integer(), sa.ForeignKey("dataset.id", ondelete="CASCADE"), index=True),
        sa.Column('uid', UUID(), server_default=text("uuid_generate_v4()"), nullable=False, index=True),
        sa.Column('type', sa.Integer(), nullable=False, index=True),
        sa.Column('status', sa.String(64), default="pending", nullable=False, index=True),
        sa.Column('name', sa.String()),
        sa.Column('description', sa.String()),
        sa.Column('input', JSONB()),
        sa.Column('output', JSONB()),
        sa.Column('location', sa.String(4096)),
        sa.Column('created_at', sa.DateTime(), default=sa.sql.func.now(), nullable=False),
    )


def downgrade():
    op.drop_table('result')
