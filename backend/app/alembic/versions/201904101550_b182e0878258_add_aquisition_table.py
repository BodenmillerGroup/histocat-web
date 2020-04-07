"""Add aquisition table

Revision ID: b182e0878258
Revises: b111fd91b7sw
Create Date: 2019-04-10 15:50:16.499555

"""
import sqlalchemy as sa
from alembic import op
from sqlalchemy.dialects.postgresql import JSONB

# revision identifiers, used by Alembic.
revision = 'b182e0878258'
down_revision = 'b111fd91b7sw'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'acquisition',
        sa.Column('id', sa.Integer(), primary_key=True, index=True),
        sa.Column('slide_id', sa.Integer(), sa.ForeignKey("slide.id", ondelete="CASCADE"), index=True),
        sa.Column('origin_id', sa.Integer(), index=True),
        sa.Column('meta', JSONB()),
        sa.Column('channels', JSONB()),
        sa.Column('location', sa.String(4096)),
        sa.Column('created_at', sa.DateTime(), default=sa.sql.func.now(), nullable=False),
    )


def downgrade():
    op.drop_table('acquisition')
