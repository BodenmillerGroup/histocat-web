"""Add dataset table

Revision ID: 56b78600703e
Revises: 59a8939464a4
Create Date: 2020-09-30 11:23:29.491360

"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy import text
from sqlalchemy.dialects.postgresql import UUID, JSONB, ARRAY


# revision identifiers, used by Alembic.
revision = '56b78600703e'
down_revision = '59a8939464a4'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'dataset',
        sa.Column('id', sa.Integer(), primary_key=True, index=True),
        sa.Column('project_id', sa.Integer(), sa.ForeignKey("project.id", ondelete="CASCADE"), index=True, nullable=False),
        sa.Column('uid', UUID(), server_default=text("uuid_generate_v4()"), index=True, nullable=False),
        sa.Column('status', sa.String(64), nullable=False, index=True),
        sa.Column('name', sa.String()),
        sa.Column('description', sa.String()),
        sa.Column('origin', sa.String()),
        sa.Column('acquisition_ids', ARRAY(sa.Integer())),
        sa.Column('channels', ARRAY(sa.String())),
        sa.Column('meta', JSONB()),
        sa.Column('location', sa.String(4096)),
        sa.Column('created_at', sa.DateTime(), default=sa.sql.func.now(), nullable=False),
    )


def downgrade():
    op.drop_table('dataset')
