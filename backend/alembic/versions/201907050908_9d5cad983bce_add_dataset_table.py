"""Add dataset table

Revision ID: 9d5cad983bce
Revises: b182e0878258
Create Date: 2019-07-05 09:08:36.605840

"""
import sqlalchemy as sa
from alembic import op
from sqlalchemy import text
from sqlalchemy.dialects.postgresql import UUID, JSONB

# revision identifiers, used by Alembic.

revision = '9d5cad983bce'
down_revision = 'b182e0878258'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'dataset',
        sa.Column('id', sa.Integer(), primary_key=True, index=True),
        sa.Column('experiment_id', sa.Integer(), sa.ForeignKey("experiment.id", ondelete="CASCADE"), index=True),
        sa.Column('user_id', sa.Integer(), sa.ForeignKey("user.id", ondelete="CASCADE"), index=True),
        sa.Column('uid', UUID(), server_default=text("uuid_generate_v4()"), nullable=False, index=True),
        sa.Column('status', sa.String(64), default="pending", nullable=False, index=True),
        sa.Column('name', sa.String()),
        sa.Column('description', sa.String()),
        sa.Column('input', JSONB()),
        sa.Column('output', JSONB()),
        sa.Column('meta', JSONB()),
        sa.Column('location', sa.String(4096)),
        sa.Column('created_at', sa.DateTime(), default=sa.sql.func.now(), nullable=False),
        sa.Column('updated_at', sa.DateTime(), default=sa.sql.func.now(), onupdate=sa.sql.func.now(), nullable=False),
    )


def downgrade():
    op.drop_table('dataset')
