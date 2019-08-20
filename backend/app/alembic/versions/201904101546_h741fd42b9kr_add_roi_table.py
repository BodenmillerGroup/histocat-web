"""Add roi table

Revision ID: h741fd42b9kr
Revises: b111fd91b7sw
Create Date: 2019-04-10 15:46:16.499555

"""
import sqlalchemy as sa
from alembic import op
from sqlalchemy.dialects.postgresql import JSONB

# revision identifiers, used by Alembic.

revision = 'h741fd42b9kr'
down_revision = 'b111fd91b7sw'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'roi',
        sa.Column('id', sa.Integer(), primary_key=True, index=True),
        sa.Column('panorama_id', sa.Integer(), sa.ForeignKey("panorama.id", ondelete="CASCADE"), index=True),
        sa.Column('metaname', sa.String(4096), index=True),
        sa.Column('original_id', sa.Integer(), index=True),
        sa.Column('meta', JSONB()),
        sa.Column('created_at', sa.DateTime(), default=sa.sql.func.now(), nullable=False),
    )


def downgrade():
    op.drop_table('roi')
