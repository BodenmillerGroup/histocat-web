"""Add roi_point table

Revision ID: c537dg62e8je
Revises: h741fd42b9kr
Create Date: 2019-04-10 15:48:16.499555

"""
import sqlalchemy as sa
from alembic import op
from sqlalchemy.dialects.postgresql import JSONB

# revision identifiers, used by Alembic.

revision = 'c537dg62e8je'
down_revision = 'h741fd42b9kr'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'roi_point',
        sa.Column('id', sa.Integer(), primary_key=True, index=True),
        sa.Column('roi_id', sa.Integer(), sa.ForeignKey("roi.id", ondelete="CASCADE"), index=True),
        sa.Column('order_number', sa.Integer()),
        sa.Column('slide_x_pos_um', sa.Float()),
        sa.Column('slide_y_pos_um', sa.Float()),
        sa.Column('panorama_pixel_x_pos', sa.Integer()),
        sa.Column('panorama_pixel_y_pos', sa.Integer()),
        sa.Column('created_at', sa.DateTime(), default=sa.sql.func.now(), nullable=False),
    )


def downgrade():
    op.drop_table('roi_point')